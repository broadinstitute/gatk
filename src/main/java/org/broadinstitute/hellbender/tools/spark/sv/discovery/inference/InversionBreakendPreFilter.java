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
 *      a) badly supported by MQ,
 *      b) mates too far away,
 *      c) covered by CPX variants
 *
 * then convert the paired mate locations into boundaries of an interval, then
 *
 * 2) filter out intervals which (see {@link OverlapBasedFilteringResult}:
 *      a) overlaps with multiple other intervals
 *      b) overlaps with no other intervals (including those that only overlaps with the "popular" ones mentioned above)
 *      c) overlaps with one other interval of the same type (INV55/55 or INV33/33 pairs)
 *
 * The resulting INV33/55 pairs are then returned via {@link OverlappingPair}.
 */
public final class InversionBreakendPreFilter {

    private enum ReasonForFilter {
        LOW_MQ, MATE_DISTANCE, COVERED_BY_CPX,
        NO_OVERLAPPER, MULTIPLE_OVERLAPPER,
        SAME_DIRECTION_OVERLAPPER // INV55/33 overlapping with INV55/33
    }

    /**
     * Holding information on two BND records
     * (INV55 and INV33, can NOT be INV33/INV33 or INV55/INV55)
     * whose mates' spanning intervals overlap.
     */
    public static final class OverlappingPair implements Serializable {
        private static final long serialVersionUID = 1L;

        final InvBreakEndContext fivePrimeBreakEnd;
        final InvBreakEndContext threePrimeBreakEnd;

        final String chr;
        final int fiveIntervalLeftBoundary;
        final int fiveIntervalRightBoundary;
        final int threeIntervalLeftBoundary;
        final int threeIntervalRightBoundary;

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

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final OverlappingPair that = (OverlappingPair) o;

            if (fiveIntervalLeftBoundary != that.fiveIntervalLeftBoundary) return false;
            if (fiveIntervalRightBoundary != that.fiveIntervalRightBoundary) return false;
            if (threeIntervalLeftBoundary != that.threeIntervalLeftBoundary) return false;
            if (threeIntervalRightBoundary != that.threeIntervalRightBoundary) return false;
            if (!fivePrimeBreakEnd.equals(that.fivePrimeBreakEnd)) return false;
            if (!threePrimeBreakEnd.equals(that.threePrimeBreakEnd)) return false;
            return chr.equals(that.chr);
        }

        @Override
        public int hashCode() {
            int result = fivePrimeBreakEnd.hashCode();
            result = 31 * result + threePrimeBreakEnd.hashCode();
            result = 31 * result + chr.hashCode();
            result = 31 * result + fiveIntervalLeftBoundary;
            result = 31 * result + fiveIntervalRightBoundary;
            result = 31 * result + threeIntervalLeftBoundary;
            result = 31 * result + threeIntervalRightBoundary;
            return result;
        }
    }

    /**
     * Main interface.
     *
     * <p>
     *     Note: this function takes lists because we know
     *     only ~500 mate pairs on one sample to be pre-processed,
     *     so paralleling with RDD may hurt performance.
     * </p>
     */
    public static List<OverlappingPair> preprocess(final List<InvBreakEndContext> inversionBreakends,
                                                   final List<VariantContext> complexVariants,
                                                   final Integer mateDistanceThreshold,
                                                   final Integer contigContigMQFilter,
                                                   final String outputPrefix,
                                                   final SAMSequenceDictionary refDict,
                                                   final Logger toolLogger) {

        final PrimitiveFilteringResult preprocessed =
                applyPrimitiveFilter(inversionBreakends, complexVariants, contigContigMQFilter, mateDistanceThreshold, refDict);

        final Comparator<InvBreakEndContext> invBreakEndContextComparator = InvBreakEndContext.makeComparator(refDict);

        final OverlapBasedFilteringResult overlapBasedFilteringResult =
                applyOverlapBasedFilter(preprocessed.picked, invBreakEndContextComparator);

        writeBedFileForFilteredVariants(outputPrefix,
                Stream.concat(preprocessed.getFilteredVariants(), overlapBasedFilteringResult.getFilteredVariants())
                        .sorted((one, two) -> one.compareTo(two, refDict)));

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
        final String bedOutput = outputPrefix + "filter.bed";
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

        PrimitiveFilteringResult(final List<InvBreakEndContext> filteredAwayDueToMQ, final List<InvBreakEndContext> filteredAwayDueToSize,
                                 final Map<InvBreakEndContext, String> coveredByComplexVariants, final SVIntervalTree<InvBreakEndContext> picked) {
            this.filteredAwayDueToMQ = filteredAwayDueToMQ;
            this.filteredAwayDueToSize = filteredAwayDueToSize;
            this.coveredByComplexVariants = coveredByComplexVariants;
            this.picked = picked;
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

    /**
     * Parse provided input VCF files,
     * split the variants by cases as designed in {@link PrimitiveFilteringResult}.
     *
     * Extract the relevant BND's
     *   turn those BND's into SVIntervals
     *   filter by evidence local assembly contigs' MQ (max of the MQ's must be >= {@code primitiveContigMQFilter})
     *   filter by size (if mate's distance is > {@code mateDistanceFilterThreshold})
     */
    private static PrimitiveFilteringResult applyPrimitiveFilter(final List<InvBreakEndContext> nonComplexVariants,
                                                                 final List<VariantContext> complexVariants,
                                                                 final int primitiveContigMQFilter,
                                                                 final int mateDistanceFilterThreshold,
                                                                 final SAMSequenceDictionary refDict) {

        final SVIntervalTree<VariantContext> complexVariantSVIntervalTree = new SVIntervalTree<>();
        complexVariants.forEach(complexVC -> {
            final SVInterval cpxInterval = new SVInterval(refDict.getSequenceIndex(complexVC.getContig()),
                    complexVC.getStart() - 1, complexVC.getEnd());
            complexVariantSVIntervalTree.put(cpxInterval, complexVC);
        });

        final List<InvBreakEndContext> filteredAwayDueToMQ = new ArrayList<>();
        final List<InvBreakEndContext> filteredAwayDueToSize = new ArrayList<>();
        final Map<InvBreakEndContext, String> coveredByComplexVariants = new HashMap<>();
        final SVIntervalTree<InvBreakEndContext> picked = new SVIntervalTree<>();
        nonComplexVariants.stream()
                .filter(SVLocalContext.BreakEndSVContext::isUpstreamMate)
                .forEach(vc -> classifyVC(vc,
                        filteredAwayDueToMQ, filteredAwayDueToSize, coveredByComplexVariants, picked,
                        mateDistanceFilterThreshold, primitiveContigMQFilter, complexVariantSVIntervalTree, refDict));

        return new PrimitiveFilteringResult(filteredAwayDueToMQ, filteredAwayDueToSize, coveredByComplexVariants, picked);
    }

    private static void classifyVC(final InvBreakEndContext invBreakEndContext,
                                   final List<InvBreakEndContext> targetContainerFilteredAwayDueToMQ,
                                   final List<InvBreakEndContext> targetContainerFilteredAwayDueToSize,
                                   final Map<InvBreakEndContext, String> targetContainerCoveredByComplexVariants,
                                   final SVIntervalTree<InvBreakEndContext> targetContainerPicked,
                                   final int mateDistanceFilterThreshold,
                                   final int primitiveContigMQFilter,
                                   final SVIntervalTree<VariantContext> complexVariantSVIntervalTree,
                                   final SAMSequenceDictionary refDict) {
        final boolean hasGoodMappings =
                invBreakEndContext.makeSureAttributeIsList(GATKSVVCFConstants.MAPPING_QUALITIES)
                        .map(Integer::valueOf).anyMatch(mq -> mq >= primitiveContigMQFilter);
        if ( ! hasGoodMappings ) {
            targetContainerFilteredAwayDueToMQ.add(invBreakEndContext);
        } else { // filter using MQ
            final SimpleInterval mateRefLoc = invBreakEndContext.getMateRefLoc();
            final SVInterval svInterval = new SVInterval(refDict.getSequenceIndex(invBreakEndContext.getContig()),
                    invBreakEndContext.getStart() - 1, mateRefLoc.getEnd());
            if (svInterval.getLength() > mateDistanceFilterThreshold)
                targetContainerFilteredAwayDueToSize.add(invBreakEndContext); // filter using size
            else {
                if (complexVariantSVIntervalTree.hasOverlapper(svInterval)) {
                    final StringBuilder stringBuilder = new StringBuilder();
                    complexVariantSVIntervalTree.overlappers(svInterval)
                            .forEachRemaining(variantContextEntry -> stringBuilder.append(variantContextEntry.getValue().getID()).append(";"));
                    targetContainerCoveredByComplexVariants.put(invBreakEndContext, stringBuilder.toString());
                } else {
                    targetContainerPicked.put(svInterval, invBreakEndContext);
                }
            }
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

    /**
     * Classify preprocessed variants as designed by {@link OverlapBasedFilteringResult}.
     */
    private static OverlapBasedFilteringResult applyOverlapBasedFilter(final SVIntervalTree<InvBreakEndContext> picked,
                                                                       final Comparator<InvBreakEndContext> invBreakEndContextComparator) {

        final Set<InvBreakEndContext> noOverlappers = new HashSet<>();
        final Map<InvBreakEndContext, List<String>> multipleOverlappers = new HashMap<>();
        final Set<Tuple2<InvBreakEndContext, InvBreakEndContext>> uniqueStrangeOverlappers = new HashSet<>();

        final Set<OverlappingPair> uniqueHetOverlappers = new HashSet<>();

        picked.forEach(e -> classifyIntervals(e, picked, noOverlappers, multipleOverlappers, uniqueStrangeOverlappers, uniqueHetOverlappers, invBreakEndContextComparator));

        // one more traverse to avoid cases that some entry's unique overlapper is records in {@code multipleOverlappers},
        // if so, reclassify that entry as having no overlapper
        postProcess(noOverlappers, multipleOverlappers, uniqueStrangeOverlappers, uniqueHetOverlappers);

        return new OverlapBasedFilteringResult(noOverlappers, multipleOverlappers, uniqueStrangeOverlappers, uniqueHetOverlappers);
    }

    private static void classifyIntervals(final SVIntervalTree.Entry<InvBreakEndContext> treeNode,
                                          final SVIntervalTree<InvBreakEndContext> pickedFromPrimitiveFiltering,
                                          final Set<InvBreakEndContext> targetContainerNoOverlappers,
                                          final Map<InvBreakEndContext, List<String>> targetContainerMultipleOverlappers,
                                          final Set<Tuple2<InvBreakEndContext, InvBreakEndContext>> targetContainerUniqueStrangeOverlappers,
                                          final Set<OverlappingPair> targetContainerUniqueHetOverlappers,
                                          final Comparator<InvBreakEndContext> invBreakEndContextComparator) {

        final SVInterval interval = treeNode.getInterval();
        final InvBreakEndContext variant = treeNode.getValue();

        final List<InvBreakEndContext> overlappers =
                Utils.stream(pickedFromPrimitiveFiltering.overlappers(interval))
                        .map(SVIntervalTree.Entry::getValue)
                        .filter(other -> !variant.equals(other)) // remove self overlap
                        .collect(Collectors.toList());

        if (overlappers.size() == 0) {
            targetContainerNoOverlappers.add(variant);
        } else if (overlappers.size() > 1) {
            targetContainerMultipleOverlappers.put(variant, overlappers.stream().map(InvBreakEndContext::getID).collect(Collectors.toList()));
        } else {
            final InvBreakEndContext uniqueOverlapper = overlappers.iterator().next();
            final int compare = invBreakEndContextComparator.compare(variant, uniqueOverlapper);

            if (variant.isType33() == uniqueOverlapper.isType33()) { // INV55/55 or INV33/33 overlapper pair
                // sorting to avoid duplicate entries in set
                if ( compare < 0 )
                    targetContainerUniqueStrangeOverlappers.add(new Tuple2<>(variant, uniqueOverlapper));
                else
                    targetContainerUniqueStrangeOverlappers.add(new Tuple2<>(uniqueOverlapper, variant));
            } else {
                if ( compare < 0 )
                    targetContainerUniqueHetOverlappers.add(new OverlappingPair(variant, uniqueOverlapper));
                else
                    targetContainerUniqueHetOverlappers.add(new OverlappingPair(uniqueOverlapper, variant));
            }
        }
    }

    private static void postProcess(final Set<InvBreakEndContext> noOverlappers,
                                    final Map<InvBreakEndContext, List<String>> multipleOverlappers,
                                    final Set<Tuple2<InvBreakEndContext, InvBreakEndContext>> uniqueStrangeOverlappers,
                                    final Set<OverlappingPair> uniqueHetOverlappers) {
        final Set<InvBreakEndContext> popularOnes = multipleOverlappers.keySet();

        for (final OverlappingPair pair : uniqueHetOverlappers) {
            if (popularOnes.contains(pair.fivePrimeBreakEnd)) {
                multipleOverlappers.get(pair.fivePrimeBreakEnd).add(pair.threePrimeBreakEnd.getID());
                noOverlappers.add(pair.threePrimeBreakEnd);
                uniqueHetOverlappers.remove(pair);
            }
            if (popularOnes.contains(pair.threePrimeBreakEnd)) {
                multipleOverlappers.get(pair.threePrimeBreakEnd).add(pair.fivePrimeBreakEnd.getID());
                noOverlappers.add(pair.fivePrimeBreakEnd);
                uniqueHetOverlappers.remove(pair);
            }
        }

        for (final Tuple2<InvBreakEndContext, InvBreakEndContext> pair : uniqueStrangeOverlappers) {
            if (popularOnes.contains(pair._1)) {
                multipleOverlappers.get(pair._1).add(pair._2.getID());
                noOverlappers.add(pair._2);
                uniqueStrangeOverlappers.remove(pair);
            }
            if (popularOnes.contains(pair._2)) {
                multipleOverlappers.get(pair._2).add(pair._1.getID());
                noOverlappers.add(pair._1);
                uniqueStrangeOverlappers.remove(pair);
            }
        }
    }

}
