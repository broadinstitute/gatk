package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.ChimericAlignment;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.GappedAlignmentSplitter;
import org.broadinstitute.hellbender.tools.spark.sv.utils.FileUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH;
import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.DiscoverVariantsFromContigAlignmentsSAMSpark.SAMFormattedContigAlignmentParser;


/**
 * This tool filters an input of SAM file containing alignments of single-ended long read
 * (be it long read sequencing, or contigs assembled from standard Illumina short reads),
 * that aims at providing an "optimal coverage" of the long read, based on an heuristic scoring scheme,
 * and saves the result as a text file, formatted as
 * (CTG_NAME, [[LIST_OF_PICKED_ALIGNMENT_INTERVAL_AS_COMPACT_STRING]]).
 */
@CommandLineProgramProperties(summary="Filters a SAM file containing long reads alignments, and outputs the filter-passing stripped down alignment information.",
        oneLineSummary="Filters a long read SAM file, and outputs essential results.",
        usageExample = "InternalFilterLongReadAlignmentsSAMSpark -I /path/to/my/dir/longReads.sam -O /path/to/my/dir/filteredAlignmentsDir",
        programGroup = StructuralVariationSparkProgramGroup.class)
@BetaFeature
public final class InternalFilterLongReadAlignmentsSAMSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(InternalFilterLongReadAlignmentsSAMSpark.class);

    @Argument(doc = "file containing non-canonical contig names (e.g chrUn_KI270588v1) in the reference, human reference assumed when omitted",
            shortName = "nonCanoTigFile",
            fullName = "nonCanonicalContigNamesFile", optional = true)
    public String nonCanonicalContigNamesFile;

    @Argument(doc = "prefix for output text file for filtered alignment intervals",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputFilePrefix;

    @Argument(doc = "whether to run old way of filtering or not",
            shortName = "OT",
            fullName = "oldFilteringToo", optional = true)
    private boolean runOldFilteringToo = false;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public ReadFilter makeReadFilter() {
        return ReadFilterLibrary.MAPPED;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaRDD<GATKRead> reads = getReads();
        final SAMFileHeader header = getHeaderForReads();

        FileUtils.writeLinesToSingleFile(
                filterByScore(reads, header, nonCanonicalContigNamesFile, localLogger)
                        .sortBy(tig -> tig.contigName, true, reads.getNumPartitions()/100) // num partition is purely guess
                        .mapToPair(contig -> new Tuple2<>(contig.contigName,
                                contig.alignmentIntervals.stream().map(AlignmentInterval::toPackedString).collect(Collectors.toList())))
                        .map(InternalFilterLongReadAlignmentsSAMSpark::formatContigInfo).collect().iterator(),
                outputFilePrefix + "_newFiltering.ai");

        if (runOldFilteringToo) {
            FileUtils.writeLinesToSingleFile(
                    oldWayOfFiltering(reads, header, localLogger).map(InternalFilterLongReadAlignmentsSAMSpark::formatContigInfo)
                    .collect().iterator(),
                    outputFilePrefix + "_oldFiltering.ai");
        }
    }

    static String formatContigInfo(final Tuple2<String, List<String>> pair) {
        return "(" + pair._1 + ",[" + pair._2 + "])";
    }

    /**
     * Delegates to {@link ChimericAlignment#parseOneContig(AlignedContig, int)}, which is currently used in SV discovery pipeline,
     * to filter out alignments and produces {@link ChimericAlignment} for variant discovery and interpretation.
     * Here it is simply appended with a collection operation that collects the alignments stored in the {@link ChimericAlignment}'s.
     */
    private static JavaPairRDD<String, List<String>> oldWayOfFiltering(final JavaRDD<GATKRead> longReads,
                                                                       final SAMFileHeader header,
                                                                       final Logger toolLogger) {

        // parse SAM records transform to AlignmentInterval format and split gapped alignment
        final JavaRDD<AlignedContig> parsedContigAlignmentsWithGapSplit
                = new SAMFormattedContigAlignmentParser(longReads, header, true, toolLogger).getAlignedContigs()
                .filter(contig -> contig.alignmentIntervals.size()>1).cache();

        // delegates to ChimericAlignment.parseOneContig()
        return
                parsedContigAlignmentsWithGapSplit
                .mapToPair(alignedContig ->
                        new Tuple2<>(alignedContig.contigName,
                                ChimericAlignment.parseOneContig(alignedContig, DEFAULT_MIN_ALIGNMENT_LENGTH).stream()
                                        .flatMap(chimericAlignment -> chimericAlignment.getAlignmentIntervals().stream())
                                        .sorted(AlignedContig.getAlignmentIntervalComparator())
                                        .collect(Collectors.toList())))
                .sortByKey()
                .mapValues(ailist -> ailist.stream().map(AlignmentInterval::toPackedString)
                                     .collect(SVUtils.arrayListCollector(ailist.size())));
    }

    /**
     * Filters an input of SAM file containing alignments of a single-ended long read that
     * aims at providing an "optimal coverage" of the long read, based on an heuristic scoring scheme
     * {@link #computeScoreOfConfiguration(List, Set, int)}.
     *
     * @param longReads     long read alignments
     * @param header        header for the long reads
     * @param toolLogger    logger for, most likely, debugging uses
     *
     * @return              contigs with alignments filtered and custom formatted as {@link AlignmentInterval}
     */
    static JavaRDD<AlignedContig> filterByScore(final JavaRDD<GATKRead> longReads,
                                                final SAMFileHeader header,
                                                final String nonCanonicalContigNamesFile,
                                                final Logger toolLogger) {

        toolLogger.info( "Processing this many raw alignments: " + longReads.count() );

        final JavaRDD<AlignedContig> parsedContigAlignments
                = new SAMFormattedContigAlignmentParser(longReads, header, false, toolLogger)
                .getAlignedContigs()
                .filter(InternalFilterLongReadAlignmentsSAMSpark::contigFilter).cache();

        toolLogger.info( "Primitive filtering based purely on MQ left these many contigs: " + parsedContigAlignments.count() );

        return filterAndSplitGappedAI(parsedContigAlignments, nonCanonicalContigNamesFile, header.getSequenceDictionary());
    }

    /**
     * Idea is to keep mapped contig that has at least two alignments over MQ 20,
     * or in the case of a single alignment, it must be MQ>20.
     */
    private static boolean contigFilter(final AlignedContig contig) {
        if (contig.alignmentIntervals.size() < 2 ) {
            return (!contig.alignmentIntervals.isEmpty()) && contig.alignmentIntervals.get(0).mapQual > 20;
        } else {
            return contig.alignmentIntervals.stream().mapToInt(ai -> ai.mapQual).filter(mq -> mq > 20).count() > 1;
        }
    }

    /**
     * Manages configuration score-based filtering,
     * then split the gapped alignments by delegating to {@link GappedAlignmentSplitter}.
     * Note that this step is essentially a flatMap operation, meaning that one contig may return 1 or multiple contigs:
     *  *) when 1 contig is yielded, it means the contig has only 1 configuration that scored the best
     *  *) when multiple contigs are yielded, it means the contig has several top-scored configurations
     * How to handle the second scenario can be treated in a separate logic unit.
     */
    @VisibleForTesting
    static JavaRDD<AlignedContig> filterAndSplitGappedAI(final JavaRDD<AlignedContig> parsedContigAlignments,
                                                         final String nonCanonicalContigNamesFile,
                                                         final SAMSequenceDictionary dictionary) {

        final Set<String> canonicalChromosomes = getCanonicalChromosomes(nonCanonicalContigNamesFile, dictionary);

        return parsedContigAlignments
                .mapToPair(alignedContig -> new Tuple2<>(alignedContig.contigName,
                        new Tuple2<>(alignedContig.contigSequence, pickBestConfigurations(alignedContig, canonicalChromosomes))))
                .flatMap(InternalFilterLongReadAlignmentsSAMSpark::reConstructContigFromPickedConfiguration);
    }

    /**
     * Pick the best configurations based on a heuristic scoring scheme implemented in
     * {@link #computeScoreOfConfiguration(List, Set, int)}.
     * @return a 2-D list, where in the case when multiple configurations are equally top-scored, all such configurations are picked up
     */
    @VisibleForTesting
    static List<List<AlignmentInterval>> pickBestConfigurations(final AlignedContig alignedContig,
                                                                final Set<String> canonicalChromosomes) {

        // group 1: get max aligner score of mappings to canonical chromosomes and speed up in case of too many mappings
        final int maxCanonicalChrAlignerScore = alignedContig.alignmentIntervals.stream()
                .filter(alignmentInterval -> canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig()))
                .mapToInt(ai -> ai.alnScore).max().orElse(0); // possible that no mapping to canonical chromosomes

        // speed up if number of alignments is too high (>10)
        // if mapped to canonical chromosomes, MQ must be >10; otherwise, must have AS higher than max canonical aligner score
        final List<AlignmentInterval> alignmentIntervals;
        if (alignedContig.alignmentIntervals.size() > 10) {
            alignmentIntervals = alignedContig.alignmentIntervals.stream()
                    .filter(alignmentInterval -> (!canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig())
                                                                           && alignmentInterval.alnScore > maxCanonicalChrAlignerScore)
                                                 || alignmentInterval.mapQual>10)
                    .collect(Collectors.toList());
        } else {
            alignmentIntervals = alignedContig.alignmentIntervals;
        }

        final int newMaxCanonicalChrAlignerScore = alignmentIntervals.stream()
                .filter(alignmentInterval -> canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig()))
                .mapToInt(ai -> ai.alnScore).max().orElse(0); // possible that no mapping to canonical chromosomes

        // group 2: generate, and score configurations
        final List<List<AlignmentInterval>> allConfigurations = Sets.powerSet(new HashSet<>(alignmentIntervals))
                .stream().map(ArrayList::new)
                // make sure within each configuration, alignments would be sorted as they would be in a corresponding AlignedContig
                .map(ls -> ls.stream().sorted(AlignedContig.getAlignmentIntervalComparator()).collect(Collectors.toList()))
                .collect(Collectors.toList());

        final List<Double> scores = allConfigurations.stream()
                .map(configuration -> computeScoreOfConfiguration(configuration, canonicalChromosomes, newMaxCanonicalChrAlignerScore))
                .collect(SVUtils.arrayListCollector(allConfigurations.size()));

        // group 3: pick the best-scored configuration(s) (if multiple configurations have equally good scores, return all of them)
        final double maxScore = scores.stream().mapToDouble(Double::doubleValue).max()
                .orElseThrow(() -> new GATKException("Cannot find best-scoring configuration on alignments of contig: " + alignedContig.contigName));

        return IntStream.range(0, allConfigurations.size())
                .filter(i -> scores.get(i) >= maxScore)
                .mapToObj(allConfigurations::get).collect(Collectors.toList());
    }

    /**
     * quick and dirty implementation for computing score of given configuration of alignments,
     * no assumption on the ordering of input alignments
     */
    @VisibleForTesting
    static double computeScoreOfConfiguration(final List<AlignmentInterval> configuration,
                                              final Set<String> canonicalChromosomes,
                                              final int maxCanonicalChrAlignerScore) {

        final double tigExplainQual = computeTigExplainQualOfOneAlignment(configuration, canonicalChromosomes, maxCanonicalChrAlignerScore);

        int redundancy = 0;
        for (int i = 0; i < configuration.size() -1 ; ++i) {
            for (int j = i + 1; j < configuration.size(); ++j) {
                final int overlap = AlignmentInterval.overlapOnContig(configuration.get(i), configuration.get(j));
                redundancy += overlap;
            }
        }

        return tigExplainQual - redundancy;
    }


    private static double computeTigExplainQualOfOneAlignment(final List<AlignmentInterval> configuration,
                                                              final Set<String> canonicalChromosomes,
                                                              final int maxCanonicalChrAlignerScore) {
        double tigExplainedQual = 0;
        for (final AlignmentInterval alignmentInterval : configuration) {
            final int len = alignmentInterval.endInAssembledContig - alignmentInterval.startInAssembledContig + 1;
            final double weight;
            if (canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig())) {
                weight = alignmentInterval.mapQual/60.0;
            } else {
                weight = Math.max(alignmentInterval.mapQual/60.0,
                                  alignmentInterval.alnScore > maxCanonicalChrAlignerScore ? 1 : 0);
            }
            tigExplainedQual += weight * len;
        }
        return tigExplainedQual;
    }

    /**
     * Depending on the number of best-scored configurations, one contig may produce multiple contigs with the same name
     * and sequence but different selected configurations.
     */
    private static Iterator<AlignedContig> reConstructContigFromPickedConfiguration(
            final Tuple2<String, Tuple2<byte[], List<List<AlignmentInterval>>>> nameSeqAndAlignmentsOfOneRead) {
        if (nameSeqAndAlignmentsOfOneRead._2._2.size() > 1) {
            return nameSeqAndAlignmentsOfOneRead._2._2.stream()
                    .map(configuration -> new AlignedContig(nameSeqAndAlignmentsOfOneRead._1, nameSeqAndAlignmentsOfOneRead._2._1,
                            splitGaps(configuration), true))
                    .sorted(sortConfigurations())
                    .collect(Collectors.toList()).iterator();
        } else {
            return Collections.singletonList(new AlignedContig(nameSeqAndAlignmentsOfOneRead._1, nameSeqAndAlignmentsOfOneRead._2._1,
                    splitGaps(nameSeqAndAlignmentsOfOneRead._2._2.get(0)), false))
                    .iterator();
        }
    }

    /**
     * when two configurations are the same, prefer the one with less alignments,
     * or less summed mismatches if still tie.
     */
    private static Comparator<AlignedContig> sortConfigurations() {
        Comparator<AlignedContig> numFirst
                = (AlignedContig x, AlignedContig y) -> Integer.compare(x.alignmentIntervals.size(), y.alignmentIntervals.size());
        Comparator<AlignedContig> mismatchSecond
                = (AlignedContig x, AlignedContig y) -> Integer.compare(x.alignmentIntervals.stream().mapToInt(ai -> ai.mismatches).sum(),
                                                                        y.alignmentIntervals.stream().mapToInt(ai -> ai.mismatches).sum());
        return numFirst.thenComparing(mismatchSecond);
    }

    private static List<AlignmentInterval> splitGaps(final List<AlignmentInterval> configuration) {
        return configuration.stream()
                .map(ai -> GappedAlignmentSplitter.split(ai, GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY,
                        SvCigarUtils.getUnclippedReadLength(ai.cigarAlong5to3DirectionOfContig)))
                .flatMap(Utils::stream).collect(Collectors.toList());
    }

    /**
     * Primary reference contigs are defined as chromosomes 1-22, X, Y, M, and defined for both GRCh38 and hg19.
     */
    static Set<String> getCanonicalChromosomes(final String nonCanonicalContigNamesFile, final SAMSequenceDictionary dictionary) {
        if (nonCanonicalContigNamesFile!= null) {

            try (final Stream<String> nonCanonical = Files.lines(Paths.get(nonCanonicalContigNamesFile))) {
                return Sets.difference(dictionary.getSequences().stream().map(SAMSequenceRecord::getSequenceName).collect(Collectors.toSet()),
                        nonCanonical.collect(Collectors.toSet()));
            } catch ( final IOException ioe ) {
                throw new UserException("Can't read nonCanonicalContigNamesFile file "+nonCanonicalContigNamesFile, ioe);
            }
        } else {
            final List<String> first22ChromosomesNum = IntStream.range(0, 23).mapToObj(String::valueOf).collect(Collectors.toList());
            final Set<String> canonicalChromosomeNames = first22ChromosomesNum.stream().map(name -> "chr" + name).collect(Collectors.toSet());
            canonicalChromosomeNames.addAll(first22ChromosomesNum);
            canonicalChromosomeNames.addAll(Arrays.asList("chrX", "chrY", "chrM", "X", "Y", "MT"));
            return canonicalChromosomeNames;
        }
    }
}
