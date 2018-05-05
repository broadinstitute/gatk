package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoverFromLocalAssemblyContigAlignmentsSpark;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import javax.annotation.Nonnull;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY;

public class AssemblyContigAlignmentsConfigPicker {

    /**
     * A filter that is used to remove contigs upfront which doesn't meet the following criteria
     * either:
     *  has only 1 mapping, with MQ strictly above this threshold
     * or:
     *  has more than 1 mapping, but only 1 mapping has MQ strictly above this threshold and it has a large gap in it.
     */
    static final int ALIGNMENT_MQ_THRESHOLD = 20;

    /**
     * A filter to boost configuration scoring implemented here:
     * if the configuration has more than 10 mappings, then
     * any mappings in such configuration with MQ
     * not strictly above this threshold is classified as bad and filtered.
     */
    static final int ALIGNMENT_MQ_THRESHOLD_FOR_SPEED_BOOST = 10;

    /**
     * Default value to be passed to
     * {@link #filterSecondaryConfigurationsByMappingQualityThreshold(List, int)}
     */
    static final int mqThreshold = 0;

    /**
     * parameters to be passed to {@link #removeNonUniqueMappings(AssemblyContigWithFineTunedAlignments, int, int)}
     * for dropping alignments that offer either low reference or read uniqueness.
     */
    public static final int ALIGNMENT_LOW_REF_UNIQUENESS_THRESHOLD = 20;
    public static final int ALIGNMENT_LOW_READ_UNIQUENESS_THRESHOLD = 10;

    /**
     * Filters an input of SAM file containing alignments of a single-ended long read that
     * aims at providing an "optimal coverage" of the assembly contig, based on an heuristic scoring scheme
     * {@link #computeScoreOfConfiguration(List, Set, int)}.
     *
     * @param assemblyAlignments    long read alignments
     * @param header                header for the long reads
     * @param scoreDiffTolerance    a tolerance where if two configurations' scores differ by less than or equal to this amount, they are considered equally good
     * @param toolLogger            logger for, most likely, debugging uses
     *
     * @return              contigs with alignments filtered and custom formatted as {@link AlignmentInterval}
     */
    public static JavaRDD<AssemblyContigWithFineTunedAlignments> createOptimalCoverageAlignmentSetsForContigs(final JavaRDD<GATKRead> assemblyAlignments,
                                                                                                              final SAMFileHeader header,
                                                                                                              final Set<String> canonicalChromosomes,
                                                                                                              final Double scoreDiffTolerance,
                                                                                                              final Logger toolLogger) {

        final JavaRDD<AlignedContig> parsedContigAlignments =
                convertRawAlignmentsToAlignedContigAndFilterByQuality(assemblyAlignments, header, toolLogger);

        final JavaPairRDD<Tuple2<String, byte[]>, List<GoodAndBadMappings>> assemblyContigWithPickedConfigurations =
                gatherBestConfigurationsForOneContig(parsedContigAlignments, canonicalChromosomes, scoreDiffTolerance);

        return assemblyContigWithPickedConfigurations
                .flatMap(AssemblyContigAlignmentsConfigPicker::reConstructContigFromPickedConfiguration)
                .map(tig -> removeNonUniqueMappings(tig, ALIGNMENT_LOW_REF_UNIQUENESS_THRESHOLD,
                        ALIGNMENT_LOW_READ_UNIQUENESS_THRESHOLD));
    }

    //==================================================================================================================

    /**
     * Parses input alignments into custom {@link AlignmentInterval} format, and
     * performs a primitive filtering implemented in
     * {@link #notDiscardForBadMQ(AlignedContig)} that
     * gets rid of contigs with no good alignments.
     */
    private static JavaRDD<AlignedContig> convertRawAlignmentsToAlignedContigAndFilterByQuality(final JavaRDD<GATKRead> assemblyAlignments,
                                                                                                final SAMFileHeader header,
                                                                                                final Logger toolLogger) {
        assemblyAlignments.cache();
        toolLogger.info( "Processing " + assemblyAlignments.count() + " raw alignments from " +
                         assemblyAlignments.map(GATKRead::getName).distinct().count() + " contigs.");

        final JavaRDD<AlignedContig> parsedContigAlignments =
                new SvDiscoverFromLocalAssemblyContigAlignmentsSpark.SAMFormattedContigAlignmentParser(assemblyAlignments, header, false)
                        .getAlignedContigs()
                        .filter(AssemblyContigAlignmentsConfigPicker::notDiscardForBadMQ).cache();
        assemblyAlignments.unpersist();
        toolLogger.info( "Filtering on MQ left " + parsedContigAlignments.count() + " contigs.");
        return parsedContigAlignments;
    }

    /**
     * Idea is to keep mapped contig that
     *  either has at least two alignments over {@link #ALIGNMENT_MQ_THRESHOLD},
     *  or in the case of a single alignment, it must be MQ > {@link #ALIGNMENT_MQ_THRESHOLD}.
     * Note that we are not simply filtering out contigs with only 1 alignment because
     * they might contain large (> 50) gaps hence should be kept.
     *
     * a point that could use improvements:
     *   the current implementation exhaustively checks the power set of all possible alignments of each assembly contig,
     *   which is computationally impossible for contigs having many-but-barely-any-good alignments, yet bringing in no value,
     *   hence this primitive filtering step to get rid of these bad assembly contigs.
     */
    @VisibleForTesting
    static boolean notDiscardForBadMQ(final AlignedContig contig) {
        if (contig.getAlignments().size() < 2 ) {
            return (!contig.getAlignments().isEmpty()) && contig.getAlignments().get(0).mapQual > ALIGNMENT_MQ_THRESHOLD;
        } else {
            int notBadMappingsCount = 0; // either more than 1 non-bad mappings, or at least 1 non-bad mapping containing large gap
            for (final AlignmentInterval alignment : contig.getAlignments()) {
                if (alignment.mapQual > ALIGNMENT_MQ_THRESHOLD) {
                    if (alignment.containsGapOfEqualOrLargerSize(GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY)) {
                        return true;// early return when a not-bad one contains a large gap
                    } else {
                        ++notBadMappingsCount;
                    }
                }
            }
            return notBadMappingsCount > 1;
        }
    }

    //==================================================================================================================

    /**
     * For each assembly contig, scores its alignment configurations and pick the best one(s).
     */
    @VisibleForTesting
    static JavaPairRDD<Tuple2<String, byte[]>, List<GoodAndBadMappings>> gatherBestConfigurationsForOneContig(final JavaRDD<AlignedContig> parsedContigAlignments,
                                                                                                              final Set<String> canonicalChromosomes,
                                                                                                              final Double scoreDiffTolerance) {

        return parsedContigAlignments
                .mapToPair(alignedContig -> new Tuple2<>(new Tuple2<>(alignedContig.getContigName(), alignedContig.getContigSequence()),
                        pickBestConfigurations(alignedContig, canonicalChromosomes, scoreDiffTolerance)))
                .mapToPair(nameSeqAndConfigurations -> new Tuple2<>(nameSeqAndConfigurations._1,
                        filterSecondaryConfigurationsByMappingQualityThreshold(nameSeqAndConfigurations._2, mqThreshold)));
    }

    /**
     * After configuration scoring and picking, the original alignments can be classified as
     * good and bad mappings:
     * good: the ones used the picked configuration
     * bad: unused alignments in the chosen configuration; these likely contain more noise than information
     *      they can be turned into string representation following the format as in {@link AlignmentInterval#toPackedString()}
     *
     * Note that a special case needs attention:
     *     if {@link #getMayBeNullGoodMappingToNonCanonicalChromosome()} returns non-null result,
     *     it is indicating an equally good--or better--non-chimeric mapping to a non-canonical chromosome exists,
     *     but to preserve the SV signal, we keep the chimeric alignments to canonical chromosomes and
     *     signal the situation to downstream units.
     */
    @VisibleForTesting
    public static final class GoodAndBadMappings {

        private final List<AlignmentInterval> goodMappings;
        private final List<AlignmentInterval> badMappings;
        private final AlignmentInterval goodMappingToNonCanonicalChromosome;

        public GoodAndBadMappings(@Nonnull final List<AlignmentInterval> goodMappings) {
            this(goodMappings, Collections.emptyList(), null);
        }

        public GoodAndBadMappings(@Nonnull final List<AlignmentInterval> goodMappings, @Nonnull final List<AlignmentInterval> badMappings,
                                  final AlignmentInterval goodMappingToNonCanonicalChr) {
            this.goodMappings = goodMappings;
            this.badMappings = badMappings;

            this.goodMappingToNonCanonicalChromosome = goodMappingToNonCanonicalChr;
        }

        public GoodAndBadMappings(@Nonnull final List<AlignmentInterval> goodMappings, @Nonnull final List<AlignmentInterval> badMappings) {
            this(goodMappings, badMappings, null);
        }

        public List<AlignmentInterval> getGoodMappings() {
            return goodMappings;
        }

        public List<AlignmentInterval> getBadMappings() {
            return badMappings;
        }

        public List<String> getBadMappingsAsCompactStrings() {
            return badMappings.stream().map(AlignmentInterval::toPackedString).collect(Collectors.toList());
        }

        public AlignmentInterval getMayBeNullGoodMappingToNonCanonicalChromosome() {
            return goodMappingToNonCanonicalChromosome;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final GoodAndBadMappings that = (GoodAndBadMappings) o;

            if (!goodMappings.equals(that.goodMappings)) return false;
            if (!badMappings.equals(that.badMappings)) return false;
            return goodMappingToNonCanonicalChromosome != null ? goodMappingToNonCanonicalChromosome.equals(that.goodMappingToNonCanonicalChromosome) : that.goodMappingToNonCanonicalChromosome == null;
        }

        @Override
        public int hashCode() {
            int result = goodMappings.hashCode();
            result = 31 * result + badMappings.hashCode();
            result = 31 * result + (goodMappingToNonCanonicalChromosome != null ? goodMappingToNonCanonicalChromosome.hashCode() : 0);
            return result;
        }

        @Override
        public String toString() {
            final StringBuilder sb = new StringBuilder("GoodAndBadMappings{");
            sb.append("goodMappings=").append(goodMappings);
            sb.append(", badMappings=").append(badMappings);
            if (goodMappingToNonCanonicalChromosome != null)
                sb.append(", goodMappingToNonCanonicalChromosome=").append(goodMappingToNonCanonicalChromosome);
            sb.append('}');
            return sb.toString();
        }
    }

    /**
     * Pick the best configurations based on a heuristic scoring scheme implemented in
     * {@link #computeScoreOfConfiguration(List, Set, int)}.
     *
     * Note that a special chanel exists, and explained in
     * {@link #getBetterNonCanonicalMapping(Set, List, int)}.
     *
     * @return a 2-D list, where in the case when multiple configurations are equally top-scored, all such configurations are picked up
     */
    @VisibleForTesting
    public static List<GoodAndBadMappings> pickBestConfigurations(final AlignedContig alignedContig,
                                                                  final Set<String> canonicalChromosomes,
                                                                  final Double scoreDiffTolerance) {
        // nothing to score if only one alignment
        if (alignedContig.getAlignments().size() == 1) {
            return Collections.singletonList(
                    new GoodAndBadMappings(Collections.singletonList(alignedContig.getAlignments().get(0)),
                                           Collections.emptyList())
            );
        }

        // step 1: get max aligner score of mappings to canonical chromosomes and speed up in case of too many mappings
        final int maxCanonicalChrAlignerScore = alignedContig.getAlignments().stream()
                .filter(alignmentInterval -> canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig()))
                .mapToInt(ai -> ai.alnScore).max().orElse(0); // possible that no mapping to canonical chromosomes

        final GoodAndBadMappings preFilteredAlignments =
                speedUpWhenTooManyMappings(alignedContig, canonicalChromosomes, maxCanonicalChrAlignerScore);
        final List<AlignmentInterval> goodMappings = preFilteredAlignments.goodMappings;
        final List<AlignmentInterval> badMappings = preFilteredAlignments.badMappings;

        final int newMaxCanonicalChrAlignerScore = goodMappings.stream()
                .filter(alignmentInterval -> canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig()))
                .mapToInt(ai -> ai.alnScore).max().orElse(0); // possible that no mapping to canonical chromosomes

        // special chanel for a special case
        final AlignmentInterval goodMappingToNonCanonicalChromosome =
                getBetterNonCanonicalMapping(canonicalChromosomes, goodMappings, newMaxCanonicalChrAlignerScore);
        if (goodMappingToNonCanonicalChromosome != null) { // take it out of consideration
            goodMappings.remove(goodMappingToNonCanonicalChromosome);
        }

        // step 2: generate, and score configurations
        return getGoodAndBadMappings(goodMappings, badMappings, goodMappingToNonCanonicalChromosome, canonicalChromosomes,
                newMaxCanonicalChrAlignerScore, scoreDiffTolerance, alignedContig.getContigName());
    }

    // speed up if number of alignments is too high (>10)
    // if mapped to canonical chromosomes, MQ must be >10; otherwise, must have AS higher than max canonical aligner score
    @VisibleForTesting
    static GoodAndBadMappings speedUpWhenTooManyMappings(final AlignedContig alignedContig,
                                                         final Set<String> canonicalChromosomes,
                                                         final int maxCanonicalChrAlignerScore) {

        final List<AlignmentInterval> goods;
        final List<AlignmentInterval> bads;
        if (alignedContig.getAlignments().size() > 10) {
            goods = new ArrayList<>();
            bads = new ArrayList<>();
            for(final AlignmentInterval alignment : alignedContig.getAlignments()) {
                final boolean isGood = (!canonicalChromosomes.contains(alignment.referenceSpan.getContig()) && alignment.alnScore > maxCanonicalChrAlignerScore)
                        || alignment.mapQual > ALIGNMENT_MQ_THRESHOLD_FOR_SPEED_BOOST;
                if (isGood)
                    goods.add(alignment);
                else
                    bads.add(alignment);
            }
        } else {
            goods = alignedContig.getAlignments();
            bads = Collections.emptyList();
        }
        return new GoodAndBadMappings(goods, bads);
    }

    private static List<GoodAndBadMappings> getGoodAndBadMappings(final List<AlignmentInterval> goodMappings,
                                                                  final List<AlignmentInterval> badMappings,
                                                                  final AlignmentInterval goodMappingToNonCanonicalChromosome,
                                                                  final Set<String> canonicalChromosomes,
                                                                  final int maxCanonicalChrAlignerScore,
                                                                  final Double scoreDiffTolerance,
                                                                  final String contigName) {
        final List<List<AlignmentInterval>> allConfigurations = Sets.powerSet(new HashSet<>(goodMappings))
                .stream().map(ArrayList::new)
                // make sure within each configuration, alignments would be sorted as they would be in a corresponding AlignedContig
                .map(ls -> ls.stream().sorted(AlignedContig.getAlignmentIntervalComparator()).collect(Collectors.toList()))
                .collect(Collectors.toList());

        final List<Double> scores = allConfigurations.stream()
                .map(configuration -> computeScoreOfConfiguration(configuration, canonicalChromosomes, maxCanonicalChrAlignerScore))
                .collect(SVUtils.arrayListCollector(allConfigurations.size()));

        // step 3: pick the best-scored configuration(s) (if multiple configurations have equally good scores, return all of them)
        final double maxScore = scores.stream().mapToDouble(Double::doubleValue).max()
                .orElseThrow(() -> new GATKException("Cannot find best-scoring configuration on alignments of contig: " + contigName));

        return IntStream.range(0, allConfigurations.size())
                .filter(i -> {
                    final Double s = scores.get(i);
                    // two configurations with would-be-same-scores can differ by a tolerance due to the sin of comparing
                    // could-be-close floating point values
                    // (see http://www.cygnus-software.com/papers/comparingfloats/Comparing%20floating%20point%20numbers.htm)
                    final Double tol = Math.max(Math.ulp(s), scoreDiffTolerance);
                    return s >= maxScore || maxScore - s <= tol;
                })
                .mapToObj(p -> {
                    final ArrayList<AlignmentInterval> copy = new ArrayList<>(goodMappings);
                    final List<AlignmentInterval> pickedAlignments = allConfigurations.get(p);
                    copy.removeAll(pickedAlignments); // remove picked, left are bad
                    copy.addAll(badMappings); // add original bad mappings
                    return new GoodAndBadMappings(pickedAlignments, copy, goodMappingToNonCanonicalChromosome);
                })
                .collect(Collectors.toList());
    }

    /**
     * quick and dirty implementation for computing score of given configuration of alignments,
     * no assumption on the ordering of input alignments
     */
    @VisibleForTesting
    static double computeScoreOfConfiguration(final List<AlignmentInterval> configuration,
                                              final Set<String> canonicalChromosomes,
                                              final int maxCanonicalChrAlignerScore) {

        final double tigExplainQual = computeTigExplainQualOfOneConfiguration(configuration, canonicalChromosomes, maxCanonicalChrAlignerScore);

        int redundancy = 0;
        for (int i = 0; i < configuration.size() -1 ; ++i) {
            for (int j = i + 1; j < configuration.size(); ++j) {
                final int overlap = AlignmentInterval.overlapOnContig(configuration.get(i), configuration.get(j));
                redundancy += overlap;
            }
        }

        return tigExplainQual - redundancy;
    }

    static final double COVERAGE_MQ_NORMALIZATION_CONST = 60.0;

    private static double computeTigExplainQualOfOneConfiguration(final List<AlignmentInterval> configuration,
                                                                  final Set<String> canonicalChromosomes,
                                                                  final int maxCanonicalChrAlignerScore) {
        double tigExplainedQual = 0;
        for (final AlignmentInterval alignmentInterval : configuration) {
            final int len = alignmentInterval.getSizeOnRead();
            final double weight;
            if (canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig())) {
                weight = alignmentInterval.mapQual/COVERAGE_MQ_NORMALIZATION_CONST;
            } else {
                weight = Math.max(alignmentInterval.mapQual/COVERAGE_MQ_NORMALIZATION_CONST,
                                  alignmentInterval.alnScore > maxCanonicalChrAlignerScore ? 1 : 0);
            }
            tigExplainedQual += weight * len;
        }
        return tigExplainedQual;
    }

    //==================================================================================================================

    /**
     * Reconstructs (possibly more than one) {@link AlignedContig} based on
     * the given best-scored configuration(s) in {@code nameSeqAndBestConfigurationsOfOneAssemblyContig}.
     *
     * @param nameSeqAndBestConfigurationsOfOneAssemblyContig the name, sequence, and picked best alignment configuration(s) of an assembly contig
     * @return The number of returned contigs will be the same as the given best-scored configurations.
     */
    @VisibleForTesting
    public static Iterator<AssemblyContigWithFineTunedAlignments> reConstructContigFromPickedConfiguration(
            final Tuple2<Tuple2<String, byte[]>, List<GoodAndBadMappings>> nameSeqAndBestConfigurationsOfOneAssemblyContig) {

        final String contigName = nameSeqAndBestConfigurationsOfOneAssemblyContig._1._1;
        final byte[] contigSeq = nameSeqAndBestConfigurationsOfOneAssemblyContig._1._2;
        final List<GoodAndBadMappings> bestConfigurations = nameSeqAndBestConfigurationsOfOneAssemblyContig._2;
        if (bestConfigurations.size() > 1) { // more than one best configuration
            return bestConfigurations.stream()
                    .map(configuration -> updateContigMappingsWithGapSplit(contigName, contigSeq, configuration, true))
                    .sorted(getConfigurationComparator())
                    .iterator();
        } else {
            return Collections.singletonList(updateContigMappingsWithGapSplit(contigName, contigSeq, bestConfigurations.get(0), false)).iterator();
        }
    }

    private static AssemblyContigWithFineTunedAlignments updateContigMappingsWithGapSplit(final String contigName, final byte[] contigSeq,
                                                                                          final GoodAndBadMappings configuration,
                                                                                          final boolean setResultContigAsAmbiguous) {
        final GoodAndBadMappings goodAndBadMappings;
        if ( configuration.goodMappings.stream().anyMatch(alignment -> alignment.containsGapOfEqualOrLargerSize(GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY)) ) {
            goodAndBadMappings = splitGaps(configuration);
        } else {
            goodAndBadMappings = configuration;
        }

        return new AssemblyContigWithFineTunedAlignments(
                new AlignedContig(contigName, contigSeq, goodAndBadMappings.goodMappings),
                goodAndBadMappings.badMappings.stream().map(AlignmentInterval::toPackedString).collect(Collectors.toList()),
                setResultContigAsAmbiguous,
                goodAndBadMappings.getMayBeNullGoodMappingToNonCanonicalChromosome());
    }

    /**
     * when two configurations are the same,
     * implement ordering that
     *      prefers the configuration with less alignments,
     *      and prefers the configuration with a lower number of summed mismatches in case of a tie
     */
    @VisibleForTesting
    static Comparator<AssemblyContigWithFineTunedAlignments> getConfigurationComparator() {
        final Comparator<AssemblyContigWithFineTunedAlignments> numFirst
                = Comparator.comparingInt(tig -> tig.getAlignments().size());
        final Comparator<AssemblyContigWithFineTunedAlignments> mismatchSecond
                = (AssemblyContigWithFineTunedAlignments x, AssemblyContigWithFineTunedAlignments y)
                -> Integer.compare(x.getAlignments().stream().mapToInt(ai -> ai.mismatches).sum(),
                                   y.getAlignments().stream().mapToInt(ai -> ai.mismatches).sum());
        return numFirst.thenComparing(mismatchSecond);
    }

    @VisibleForTesting
    static GoodAndBadMappings splitGaps(final GoodAndBadMappings configuration) {

        final List<AlignmentInterval> originalBadMappings = configuration.badMappings;
        final List<AlignmentInterval> scan = configuration.goodMappings;

        // 1st pass, split gapped alignments when available, and all defaults to good
        final List<Tuple2<Boolean, Iterable<AlignmentInterval>>> alignmentSplitChildren =
                scan.stream()
                        .map(alignment -> {
                            final Iterable<AlignmentInterval> split;
                            if (alignment.containsGapOfEqualOrLargerSize(GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY)) {
                                split = ContigAlignmentsModifier.splitGappedAlignment(alignment, GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY,
                                        SvCigarUtils.getUnclippedReadLength(alignment.cigarAlong5to3DirectionOfContig));
                            } else {
                                split = Collections.singletonList(alignment);
                            }
                            return new Tuple2<>(true, split);
                        }).collect(Collectors.toList());

        // 2nd pass make a choice between gapped and overlapping alignment (alignments that are not favored has its "2nd" set to false)
        final int count = scan.size();
        for (int i = 0; i < count; ++i) {
            final AlignmentInterval alignment = scan.get(i);
            final Tuple2<Boolean, Iterable<AlignmentInterval>> split = alignmentSplitChildren.get(i);
            if ( split._1 && Iterables.size(split._2) != 1 ) { // the split could be already marked bad (i.e. to be filtered out), or could contain no large gaps (i.e. should not check it)
                for (int j = 0; j < count; ++j) {
                    final AlignmentInterval other = scan.get(j);
                    if (j == i || AlignmentInterval.overlapOnContig(alignment, other) == 0)
                        continue;

                    if ( Utils.stream(split._2).anyMatch(other::containsOnRead) ) {
                        if ( gappedAlignmentOffersBetterCoverage(alignment, other) ) {
                            final Iterable<AlignmentInterval> copy = alignmentSplitChildren.get(j)._2;
                            alignmentSplitChildren.set(j, new Tuple2<>(false, copy));
                        } else {
                            final Iterable<AlignmentInterval> copy = alignmentSplitChildren.get(i)._2;
                            alignmentSplitChildren.set(i, new Tuple2<>(false, copy));
                        }
                    }
                }
            }
        }

        final List<AlignmentInterval> good = new ArrayList<>();
        final List<AlignmentInterval> bad = new ArrayList<>(originalBadMappings);
        for (final Tuple2<Boolean, Iterable<AlignmentInterval>> pair : alignmentSplitChildren) {
            if (pair._1) {
                good.addAll( Lists.newArrayList(pair._2) );
            } else {
                bad.addAll( Lists.newArrayList(pair._2) );
            }
        }
        return new GoodAndBadMappings(good, bad, configuration.goodMappingToNonCanonicalChromosome);
    }

    @VisibleForTesting
    static boolean gappedAlignmentOffersBetterCoverage(final AlignmentInterval gapped,
                                                       final AlignmentInterval overlappingNonGapped) {
        final int diff = gapped.getSizeOnRead() - overlappingNonGapped.getSizeOnRead();
        if ( diff == 0) {
            return gapped.alnScore > overlappingNonGapped.alnScore;
        } else {
            return diff > 0;
        }
    }

    //==================================================================================================================

    /**
     * For contigs with more than 1 best-scored configurations as determined by
     * {@link #pickBestConfigurations(AlignedContig, Set, Double)},
     * save the contigs that has one and only one configuration that
     * has all mapping quality strictly above the specified {@code threshold}.
     * Example:
     *  if a contig has two equal scored configurations with MQ's {10, 60, 60}, and {60, 60},
     *  this function will favor/pick the {60, 60} configuration if the {@code threshold} is 10,
     *  hence remove the ambiguity;
     *  on the other hand if the {@code threshold} is passed in as any value below 10, say 0,
     *  then this function returns all the original {@code differentRepresentationsForOneContig}
     */
    @VisibleForTesting
    static List<GoodAndBadMappings> filterSecondaryConfigurationsByMappingQualityThreshold(
            final List<GoodAndBadMappings> differentRepresentationsForOneContig,
            final int threshold) {

        if ( differentRepresentationsForOneContig.size() == 1) {
            return differentRepresentationsForOneContig;
        } else {
            final List<GoodAndBadMappings> collect = Utils.stream(differentRepresentationsForOneContig)
                    .filter(rep -> rep.goodMappings.stream().mapToInt(ai -> ai.mapQual).min().orElse(threshold) > threshold)
                    .collect(Collectors.toList());
            if ( collect.size()!=1 ) {
                return differentRepresentationsForOneContig;
            } else {
                return collect;
            }
        }
    }

    /**
     * There are locations on the non-canonical chromosomes of the HG38 reference that are similar to a location
     * on the canonical chromosomes, except that it has rearranged (or deleted/duplicated) parts of the canonical chromosomes.
     * In other words, the non-canonical chromosome captures an SV--relative to the canonical chromosomes--of relatively high frequency.
     *
     * The sample under analysis could have the allele of this non-canonical version,
     * and be marked as having an SV on the corresponding location on the canonical chromosome.
     *
     * An assembly contig from this sample may have two equally well scored alignment configurations, where
     *  one configuration has split alignments to canonical chromosomes, hence indicating the SV, whereas
     *  the other configuration has a single, often very good (or even better) alignment to a non-canonical chromosome.
     * We send down the chimeric alignment configuration for inference, but notes down that a non-canonical chromosome
     * in the reference input could have already captured the SV on this sample.
     *
     * @return  {@code null} if the non-canonical chromosome mapping doesn't offer a better score,
     *          otherwise the non-canonical chromosome mapping
     */
    @VisibleForTesting
    static AlignmentInterval getBetterNonCanonicalMapping(final Set<String> canonicalChromosomes,
                                                          final List<AlignmentInterval> goodMappings,
                                                          final int maxCanonicalChrAlignerScore) {
        final List<AlignmentInterval> canonicalMappings = new ArrayList<>(goodMappings.size());
        final List<AlignmentInterval> nonCanonicalMapping = new ArrayList<>();
        for (final AlignmentInterval alignment : goodMappings) {
            if ( canonicalChromosomes.contains(alignment.referenceSpan.getContig()) )
                canonicalMappings.add(alignment);
            else
                nonCanonicalMapping.add(alignment);
        }
        if ( nonCanonicalMapping.size() == 1 &&
                (canonicalMappings.size() > 1 || canonicalMappings.get(0).containsGapOfEqualOrLargerSize(GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY)) ) {
            final double canonicalScore = computeScoreOfConfiguration(canonicalMappings, canonicalChromosomes, maxCanonicalChrAlignerScore);
            final double nonCanonicalScore = computeScoreOfConfiguration(nonCanonicalMapping, canonicalChromosomes, maxCanonicalChrAlignerScore);
            return ( canonicalScore > nonCanonicalScore ) ? null : nonCanonicalMapping.get(0);
        } else {
            return null;
        }
    }

    //==================================================================================================================

    /**
     * Process provided {@code originalConfiguration} of an assembly contig and split between good and bad alignments.
     *
     * <p>
     *     What is considered good and bad?
     *     For a particular mapping/alignment, it may offer low uniqueness in two sense:
     *     <ul>
     *         <li>
     *             low REFERENCE UNIQUENESS: meaning the sequence being mapped match multiple locations on the reference;
     *         </li>
     *         <li>
     *             low READ UNIQUENESS: with only a very short part of the read uniquely explained by this particular alignment;
     *         </li>
     *     </ul>
     *     Good alignments offer both high reference uniqueness and read uniqueness, as judged by the requested
     *     {@code mapQThresholdInclusive} and {@code uniqReadLenInclusive}
     *     (yes we are doing a hard-filtering but more advanced model is not our priority right now 2017-11-20).
     * </p>
     *
     * Note that "original" is meant to be possibly different from the returned configuration,
     * but DOES NOT mean the alignments of the contig as given by the aligner,
     * i.e. the configuration should be one of the best given by
     * {@link #pickBestConfigurations(AlignedContig, Set, Double)}.
     */
    public static AssemblyContigWithFineTunedAlignments removeNonUniqueMappings(final AssemblyContigWithFineTunedAlignments tig,
                                                                                final int mapQThresholdInclusive,
                                                                                final int uniqReadLenInclusive) {
        final List<AlignmentInterval> inputAlignments = tig.getAlignments();
        if (inputAlignments.size() <= 2) return tig;

        // two pass, each focusing on removing the alignments of a contig that offers low uniqueness in one sense:

        // first pass is for removing alignments with low REFERENCE UNIQUENESS, using low mapping quality as the criterion
        final List<AlignmentInterval> selectedAlignments = new ArrayList<>(inputAlignments.size());
        final List<AlignmentInterval> lowUniquenessMappings = new ArrayList<>(inputAlignments.size());

        for (final AlignmentInterval alignment : inputAlignments) {
            if (alignment.mapQual >= mapQThresholdInclusive)
                selectedAlignments.add(alignment);
            else
                lowUniquenessMappings.add(alignment);
        }

        // second pass, the slower one, is to remove alignments offering low READ UNIQUENESS,
        // i.e. with only a very short part of the read being uniquely explained by this particular alignment;
        // the steps are:
        //      search bi-directionally until cannot find overlap any more,
        //      subtract the overlap from the distance covered on the contig by the alignment.
        //      This gives unique read region it explains.
        //      If this unique read region is "short": shorter than {@code uniqReadLenInclusive}), drop it.

        // each alignment has an entry of a tuple2, one for max overlap maxFront, one for max overlap maxRear,
        // max overlap maxFront is a tuple2 registering the index and overlap bases count
        final Map<AlignmentInterval, Tuple2<Integer, Integer>> maxOverlapMap = getMaxOverlapPairs(selectedAlignments);
        for(Iterator<AlignmentInterval> iterator = selectedAlignments.iterator(); iterator.hasNext();) {
            final AlignmentInterval alignment = iterator.next();

            final Tuple2<Integer, Integer> maxOverlapFrontAndRear = maxOverlapMap.get(alignment);
            final int maxOverlapFront = Math.max(0, maxOverlapFrontAndRear._1);
            final int maxOverlapRear = Math.max(0, maxOverlapFrontAndRear._2);

            // theoretically this could be negative for an alignment whose maxFront and maxRear sums together bigger than the read span
            // but earlier configuration scoring would make this impossible because such alignments should be filtered out already
            // considering that it brings more penalty than value, i.e. read bases explained (even if the MQ is 60),
            // but even if it is kept, a negative value won't hurt unless a stupid threshold value is passed in
            final int uniqReadSpan = alignment.endInAssembledContig - alignment.startInAssembledContig + 1
                    - maxOverlapFront - maxOverlapRear;
            if (uniqReadSpan < uniqReadLenInclusive) {
                lowUniquenessMappings.add(alignment);
                iterator.remove();
            }
        }

        final AlignedContig updatedTig = new AlignedContig(tig.getContigName(), tig.getContigSequence(), selectedAlignments);
        return new AssemblyContigWithFineTunedAlignments(
                updatedTig,
                lowUniquenessMappings.stream().map(AlignmentInterval::toPackedString).collect(Collectors.toList()),
                tig.hasEquallyGoodAlnConfigurations(),
                tig.getSAtagForGoodMappingToNonCanonicalChromosome());
    }

    /**
     * Each alignment in a specific configuration has an entry,
     * pointing to the alignments that comes before and after it,
     * that overlaps maximally (i.e. no other front or rear alignments have more overlaps)
     * with the current alignment.
     */
    private static final class TempMaxOverlapInfo {
        final Tuple2<Integer, Integer> maxFront; // 1st holds index pointing to another alignment before this, 2nd holds the count of overlapping bases
        final Tuple2<Integer, Integer> maxRear;  // same intention as above, but for alignments after this

        TempMaxOverlapInfo() {
            maxFront = new Tuple2<>(-1, -1);
            maxRear = new Tuple2<>(-1 ,-1);
        }

        TempMaxOverlapInfo(final Tuple2<Integer, Integer> maxFront, final Tuple2<Integer, Integer> maxRear) {
            this.maxFront = maxFront;
            this.maxRear = maxRear;
        }
    }

    /**
     * Extract the max overlap information, front and back, for each alignment in {@code configuration}.
     * For each alignment, the corresponding tuple2 has the max (front, rear) overlap base counts.
     */
    private static Map<AlignmentInterval, Tuple2<Integer, Integer>> getMaxOverlapPairs(final List<AlignmentInterval> configuration) {

        final List<TempMaxOverlapInfo> intermediateResult =
                new ArrayList<>(Collections.nCopies(configuration.size(), new TempMaxOverlapInfo()));

        // We iterate through all alignments except the last one
        // For the last alignment, which naturally doesn't have any maxRear,
        //     the following implementation sets its maxFront during the iteration, if available at all (it may overlap with nothing)
        for(int i = 0; i < configuration.size() - 1; ++i) {

            final AlignmentInterval cur = configuration.get(i);
            // For the i-th alignment, we only look at alignments after it (note j starts from i+1) and find max overlap
            int maxOverlapRearBases = -1;
            int maxOverlapRearIndex = -1;
            for (int j = i + 1; j < configuration.size(); ++j) { // note j > i
                final int overlap = AlignmentInterval.overlapOnContig(cur, configuration.get(j));
                if (overlap > maxOverlapRearBases) {
                    maxOverlapRearBases = overlap;
                    maxOverlapRearIndex = j;
                } else { // following ones, as guaranteed by the ordering of alignments in the contig, cannot overlap
                    break;
                }
            }

            if (maxOverlapRearBases > 0){
                // for current alignment (i-th), set its max_overlap_rear, which would not change in later iterations and copy old max_overlap_front
                final Tuple2<Integer, Integer> maxRear = new Tuple2<>(maxOverlapRearIndex, maxOverlapRearBases);
                final Tuple2<Integer, Integer> maxFrontToCopy = intermediateResult.get(i).maxFront;
                intermediateResult.set(i, new TempMaxOverlapInfo(maxFrontToCopy, maxRear));

                // then conditionally set the max_overlap_front of the
                // maxOverlapRearIndex-th alignment
                // that maximally overlaps with the current, i.e. i-th, alignment
                final TempMaxOverlapInfo oldValue = intermediateResult.get(maxOverlapRearIndex);// maxOverlapRearIndex cannot be -1 here
                if (oldValue.maxFront._2 < maxOverlapRearBases)
                    intermediateResult.set(maxOverlapRearIndex, new TempMaxOverlapInfo(new Tuple2<>(i, maxOverlapRearBases), oldValue.maxRear));
            }
        }

        final Map<AlignmentInterval, Tuple2<Integer, Integer>> maxOverlapMap = new HashMap<>(configuration.size());
        for (int i = 0; i < configuration.size(); ++i) {
            maxOverlapMap.put(configuration.get(i),
                    new Tuple2<>(intermediateResult.get(i).maxFront._2, intermediateResult.get(i).maxRear._2));
        }

        return maxOverlapMap;
    }
}
