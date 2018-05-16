package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVLocalContext;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import scala.Tuple2;

import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.SVLocalContext.InvBreakEndContext;

/**
 * Strategy is to perform pre-filtering upfront because they are likely to be artifact or different types:
 *
 * 1) first filter out BND records that are (see {@link PrimitiveFilteringResult}):
 *      a) inter-chromosome, no SS, 2nd in mate (redundant),
 *      b) badly supported by MQ,
 *      c) mates too far away,
 *      d) covered by CPX variants
 *
 * then convert the paired mate locations into boundaries of an interval, then
 *
 * 2) filter out intervals that are (see {@link OverlapBasedFilteringResult}:
 *      a) overlaps with multiple other intervals
 *      b) overlaps with no other intervals (including those that only overlaps with the "popular" ones mentioned above)
 *      c) overlaps with one other interval of the same type (INV55/55 or INV33/33 pairs)
 *
 * The resulting INV33/55 pairs are then returned.
 */
final class InversionBreakendPreFilter {

    private enum ReasonForFilter {
        LOW_MQ, MATE_DISTANCE, COVERED_BY_CPX,
        NO_OVERLAPPER, MULTIPLE_OVERLAPPER,
        SAME_DIRECTION_OVERLAPPER // INV55/33 overlapping with INV55/33
    }

    public static final class OverlappingPair implements Serializable {
        private static final long serialVersionUID = 1L;

        protected final InvBreakEndContext fivePrimeBreakEnd;
        protected final InvBreakEndContext threePrimeBreakEnd;

        protected final String chr;
        protected final int fiveIntervalLeftBoundary;
        protected final int fiveIntervalRightBoundary;
        protected final int threeIntervalLeftBoundary;
        protected final int threeIntervalRightBoundary;

        OverlappingPair(final SVLocalContext.InvBreakEndContext first, final SVLocalContext.InvBreakEndContext second) {

            if (first.isType33()) {
                fivePrimeBreakEnd = second;
                threePrimeBreakEnd = first;
            } else {
                fivePrimeBreakEnd = first;
                threePrimeBreakEnd = second;
            }
            chr = fivePrimeBreakEnd.getContig();
            fiveIntervalLeftBoundary = fivePrimeBreakEnd.getStart();
            fiveIntervalRightBoundary = fivePrimeBreakEnd.getMateRefLoc().getEnd();
            threeIntervalLeftBoundary = threePrimeBreakEnd.getStart();
            threeIntervalRightBoundary = threePrimeBreakEnd.getMateRefLoc().getEnd();
        }
    }

    // note: this function takes lists because we know there's going to be only ~500 mate pairs on a particular sample,
    // so paralleling with RDD may hurt performance
    static List<OverlappingPair> preprocess(final List<InvBreakEndContext> inversionBreakends,
                                            final List<VariantContext> complexVariants,
                                            final Integer mateDistanceThreshold,
                                            final Integer contigContigMQFilter,
                                            final String outputPrefix,
                                            final SAMSequenceDictionary refDict,
                                            final Logger toolLogger) {

        final PrimitiveFilteringResult preprocessed;

        final Comparator<InvBreakEndContext> invBreakEndContextComparator = InvBreakEndContext.makeComparator(refDict);

        final OverlapBasedFilteringResult overlapBasedFilteringResult;

        return overlapBasedFilteringResult.uniqueHetOverlappers.stream()
                .sorted((pair1, pair2) -> {
                    int compare = invBreakEndContextComparator.compare(pair1.fivePrimeBreakEnd, pair2.fivePrimeBreakEnd);
                    if ( 0 == compare ) {
                        return invBreakEndContextComparator.compare(pair1.threePrimeBreakEnd, pair2.threePrimeBreakEnd);
                    } else
                        return compare;
                }).collect(Collectors.toList());
    }

    private static final class AnnotatedFilteredInterval {
        private final SimpleInterval interval;   // interval
        private final String annotation;

        private AnnotatedFilteredInterval(final SimpleInterval interval, final String annotation) {
            this.interval = interval;
            this.annotation = annotation;
        }

        int compareTo(final AnnotatedFilteredInterval other, final SAMSequenceDictionary refDict) {
            return IntervalUtils.compareLocatables(interval, other.interval, refDict);
        }

        private String toBedLine() {
            return interval.getContig() + "\t" + interval.getStart() + "\t" + interval.getEnd() + "\t" + annotation;
        }
    }

    private static void writeBedFileForFilteredVariants(final String outputPrefix,
                                                        final Stream<AnnotatedFilteredInterval> annotatedMateIntervals) {
        final String bedOutput = outputPrefix + ".bed";
        try {
            Files.write(IOUtils.getPath(bedOutput),
                    annotatedMateIntervals
                            .map(interval -> (CharSequence)interval.toBedLine())
                            .collect(Collectors.toList()));
        } catch (final IOException ex) {
            throw new GATKException("Cannot write BED file for filtered INV breakend pairs " + bedOutput);
        }
    }

    //==================================================================================================================

    /**
     * Note that this class holds the 1st mate only,
     * since we currently don't have records with more than one mate,
     * hence the only one mate is enough.
     */
    private static final class PrimitiveFilteringResult {
        private final List<InvBreakEndContext> filteredAwayDueToMQ;             // variants whose evidence contigs' MQ are below a specified value
        private final List<InvBreakEndContext> filteredAwayDueToSize;           // variants whose distance between mates are over a specified value
        private final Map<InvBreakEndContext, String> coveredByComplexVariants; // variants whose corresponding interval is overlapping a detected Cpx variant

        private final SVIntervalTree<InvBreakEndContext> picked;                // variants picked for downstream analysis

        private PrimitiveFilteringResult(final List<InvBreakEndContext> filteredAwayDueToMQ,
                                         final List<InvBreakEndContext> filteredAwayDueToSize,
                                         final SVIntervalTree<InvBreakEndContext> picked) {
            this.filteredAwayDueToMQ = filteredAwayDueToMQ;
            this.filteredAwayDueToSize = filteredAwayDueToSize;
            this.picked = picked;
            coveredByComplexVariants = new HashMap<>();
        }

        private void moveFromPickedIfCoveredByComplex(final VariantContext complexVC, final SAMSequenceDictionary refDict) {
            final SVInterval cpxInterval = new SVInterval(refDict.getSequenceIndex(complexVC.getContig()),
                    complexVC.getStart() - 1, complexVC.getEnd());
            if ( picked.hasOverlapper(cpxInterval) ) {
                picked.overlappers(cpxInterval)
                        .forEachRemaining(bnd -> {
                            picked.remove(bnd.getInterval());
                            coveredByComplexVariants.put(bnd.getValue(), complexVC.getID());
                        });
            }
        }

        private Stream<AnnotatedFilteredInterval> getFilteredVariants() {

            Stream<AnnotatedFilteredInterval> mq = filteredAwayDueToMQ.stream().map(vc ->
                    new AnnotatedFilteredInterval(new SimpleInterval(vc.getContig(), vc.getStart(), vc.getMateRefLoc().getEnd()),
                            vc.getID().replaceAll("_1$", "") + ";" + ReasonForFilter.LOW_MQ.name())
            );
            Stream<AnnotatedFilteredInterval> dist = filteredAwayDueToSize.stream().map(vc ->
                    new AnnotatedFilteredInterval(new SimpleInterval(vc.getContig(), vc.getStart(), vc.getMateRefLoc().getEnd()),
                            vc.getID().replaceAll("_1$", "") + ";" + ReasonForFilter.MATE_DISTANCE.name())
            );
            Stream<AnnotatedFilteredInterval> covered = coveredByComplexVariants.entrySet().stream().map(kV -> {
                        InvBreakEndContext vc = kV.getKey();
                        return new AnnotatedFilteredInterval(new SimpleInterval(vc.getContig(), vc.getStart(), vc.getMateRefLoc().getEnd()),
                                vc.getID().replaceAll("_1$", "") + ";" + ReasonForFilter.COVERED_BY_CPX + ":" + kV.getValue());
            });

            return Stream.of(mq, dist, covered).flatMap(Function.identity());
        }
    }


    //==================================================================================================================

    /**
     * Further cuts down the cases to look at after applying filters in {@link PrimitiveFilteringResult},
     * that is, "filter" out breakends that
     *   1) is not paired, i.e. the interval bounded by the suggested novel adjacency does not overlap with any other intra-chromosome SS BND s;
     *   2) is paired with multiple intra-chromosome SS BND s;
     *   3) is paired with one intra-chromosome SS BND, but a INV55/INV55 or INV33/INV33 pair
     */
    private static final class OverlapBasedFilteringResult {
        private final Set<InvBreakEndContext> noOverlappers;                                        // variants whose has no overlapping INV BND suspects, or whose unique overlapper is one of {@code multipleOverlappers}
        private final Map<InvBreakEndContext, List<String>> multipleOverlappers;                    // variants who overlaps multiple other records (note: not using map because VC.equals() is hard)
        private final Set<Tuple2<InvBreakEndContext, InvBreakEndContext>> uniqueStrangeOverlappers; // pairs of variants (INV55/55, INV33/33) that overlap each other uniquely

        private final Set<OverlappingPair> uniqueHetOverlappers;                                    // pairs of variants (INV55/33, INV33/55) that overlap each other uniquely, pairs that will be analyzed

        OverlapBasedFilteringResult(final Set<InvBreakEndContext> noOverlappers,
                                    final Map<InvBreakEndContext, List<String>> multipleOverlappers,
                                    final Set<Tuple2<InvBreakEndContext, InvBreakEndContext>> uniqueStrangeOverlappers,
                                    final Set<OverlappingPair> uniqueHetOverlappers) {
            this.noOverlappers = noOverlappers;
            this.multipleOverlappers = multipleOverlappers;
            this.uniqueStrangeOverlappers = uniqueStrangeOverlappers;
            this.uniqueHetOverlappers = uniqueHetOverlappers;
        }

        private Stream<AnnotatedFilteredInterval> getFilteredVariants() {

            Stream<AnnotatedFilteredInterval> lonely = noOverlappers.stream().map(vc ->
                    new AnnotatedFilteredInterval(new SimpleInterval(vc.getContig(), vc.getStart(), vc.getMateRefLoc().getEnd()),
                            vc.getID().replaceAll("_1$", "") + ";" + ReasonForFilter.NO_OVERLAPPER.name())
            );
            Stream<AnnotatedFilteredInterval> tooPopular = multipleOverlappers.entrySet().stream().map(kV -> {
                InvBreakEndContext vc = kV.getKey();
                return new AnnotatedFilteredInterval(new SimpleInterval(vc.getContig(), vc.getStart(), vc.getMateRefLoc().getEnd()),
                        vc.getID().replaceAll("_1$", "") + ";" + ReasonForFilter.MULTIPLE_OVERLAPPER + ":" + String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, kV.getValue()));
            });

            uniqueStrangeOverlappers.stream().flatMap(pair -> {
                SimpleInterval firstInterval = new SimpleInterval(pair._1.getContig(), pair._1.getStart(), pair._1.getMateRefLoc().getEnd());
                String firstAnnotation = pair._1.getID().replaceAll("_1$", "") + ";" + ReasonForFilter.MULTIPLE_OVERLAPPER + ":" + pair._2.getID().replaceAll("_1$", "");
                SimpleInterval secondInterval = new SimpleInterval(pair._2.getContig(), pair._2.getStart(), pair._2.getMateRefLoc().getEnd());
                String secondAnnotation = pair._2.getID().replaceAll("_1$", "") + ";" + ReasonForFilter.MULTIPLE_OVERLAPPER + ":" + pair._1.getID().replaceAll("_1$", "");
                return Stream.of(new AnnotatedFilteredInterval(firstInterval, firstAnnotation),
                                 new AnnotatedFilteredInterval(secondInterval, secondAnnotation));
            });

            return Stream.of(lonely, tooPopular).flatMap(Function.identity());
        }
    }

}
