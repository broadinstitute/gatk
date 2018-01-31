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
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.GappedAlignmentSplitter;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoverFromLocalAssemblyContigAlignmentsSpark;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.ChimericAlignment;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.*;


/**
 * (Internal) Examines alignments of chimeric contigs, attempting to produce an optimal tiling
 *
 * <p>This tool is used in development and should not be of interest to most researchers.  It is a prototype
 * of one aspect of structural variant calling from chimeric contigs produced by local assemblies.</p>
 * <p>It takes a SAM/BAM/CRAM containing the alignments of assembled contigs and filters them with the
 * aim of providing "optimal coverage" of the contig, based on an heuristic scoring scheme.</p>
 * <p>It saves the result as a text file, formatted as:</p>
 * <pre>(CONTIG_NAME, [[LIST_OF_PICKED_ALIGNMENT_INTERVALS_AS_A_COMPACT_STRING]])</pre>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>An input file of contigs from local assemblies aligned to reference.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A text file describing the selected alignments.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk FilterLongReadAlignmentsSAMSpark \
 *     -I assemblies.sam \
 *     -O selected_alignments.txt
 * </pre>
 * <p>This tool can be run without explicitly specifying Spark options. That is to say, the given example command
 * without Spark options will run locally. See
 * <a href ="https://software.broadinstitute.org/gatk/documentation/article?id=10060">Tutorial#10060</a>
 * for an example of how to set up and run a Spark tool on a cloud Spark cluster.</p>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "(Internal) Examines alignments of chimeric contigs, attempting to produce an optimal tiling",
        summary =
        "This tool is used in development and should not be of interest to most researchers.  It is a prototype" +
        " of one aspect of structural variant calling from chimeric contigs produced by local assemblies." +
        " It takes a SAM/BAM/CRAM containing the alignments of assembled contigs and filters them with the" +
        " aim of providing \"optimal coverage\" of the contig, based on an heuristic scoring scheme.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public final class FilterLongReadAlignmentsSAMSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private final Logger localLogger = LogManager.getLogger(FilterLongReadAlignmentsSAMSpark.class);

    @Advanced
    @Argument(doc = "Maximum difference between alignment configuration scores for which two configurations will be considered equally likely",
            shortName = "cst", fullName = "config-score-diff-tolerance", optional = true)
    public final Double configScoreDiffTolerance = 0.0;

    @Argument(doc = "file containing non-canonical contig names (e.g chrUn_KI270588v1) in the reference, human reference assumed when omitted",
            shortName = "alt-tigs", fullName = "non-canonical-contig-names-file", optional = true)
    public String nonCanonicalContigNamesFile;

    @Argument(doc = "prefix for output text file for filtered alignment intervals",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputFilePrefix;

    @Argument(doc = "whether to run old way of filtering or not",
            shortName = "OT", fullName = "old-filtering-too", optional = true)
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

        try {
            final JavaRDD<GATKRead> reads = getReads();
            final SAMFileHeader header = getHeaderForReads();

            Files.write(Paths.get(outputFilePrefix + "_newFiltering.ai"),
                    () -> filterByScore(reads, header, nonCanonicalContigNamesFile, configScoreDiffTolerance, localLogger)
                            .sortBy(tig -> tig.contigName, true, reads.getNumPartitions() / 100) // num partition is purely guess
                            .mapToPair(contig -> new Tuple2<>(contig.contigName,
                                    contig.alignmentIntervals.stream().map(AlignmentInterval::toPackedString).collect(Collectors.toList())))
                            .map(AlignedContig::formatContigInfo)
                            .map(s -> (CharSequence)s).collect().iterator()
            );

            if (runOldFilteringToo) {
                Files.write(Paths.get(outputFilePrefix + "_oldFiltering.ai"),
                        () -> oldWayOfFiltering(reads, header)
                                .map(AlignedContig::formatContigInfo)
                                .map(s -> (CharSequence)s).collect().iterator()
                );
            }
        } catch (final IOException ioe) {
            throw new UserException.CouldNotCreateOutputFile("Could not save filtering results to file", ioe);
        }
    }

    /**
     * Delegates to {@link ChimericAlignment#parseOneContig(AlignedContig, SAMSequenceDictionary, boolean, int, int, boolean)}, which is currently used in SV discovery pipeline,
     * to filter out alignments and produces {@link ChimericAlignment} for variant discovery and interpretation.
     * Here it is simply appended with a collection operation that collects the alignments stored in the {@link ChimericAlignment}'s.
     */
    private static JavaPairRDD<String, List<String>> oldWayOfFiltering(final JavaRDD<GATKRead> longReads,
                                                                       final SAMFileHeader header) {

        // parse SAM records transform to AlignmentInterval format and split gapped alignment
        final JavaRDD<AlignedContig> parsedContigAlignmentsWithGapSplit =
                new SvDiscoverFromLocalAssemblyContigAlignmentsSpark.SAMFormattedContigAlignmentParser(longReads, header, true)
                        .getAlignedContigs()
                        .filter(contig -> contig.alignmentIntervals.size()>1).cache();

        final SAMSequenceDictionary referenceDictionary = header.getSequenceDictionary();
        // delegates to ChimericAlignment.parseOneContig()
        return
                parsedContigAlignmentsWithGapSplit
                .mapToPair(alignedContig ->
                        new Tuple2<>(alignedContig.contigName,
                                ChimericAlignment.parseOneContig(alignedContig, referenceDictionary,
                                        true, DEFAULT_MIN_ALIGNMENT_LENGTH,
                                        CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD, true).stream()
                                        .flatMap(chimericAlignment -> chimericAlignment.getAlignmentIntervals().stream())
                                        .sorted(AlignedContig.getAlignmentIntervalComparator())
                                        .collect(Collectors.toList())))
                .sortByKey()
                .mapValues(aiList -> aiList.stream().map(AlignmentInterval::toPackedString)
                                     .collect(SVUtils.arrayListCollector(aiList.size())));
    }

    /**
     * Filters an input of SAM file containing alignments of a single-ended long read that
     * aims at providing an "optimal coverage" of the long read, based on an heuristic scoring scheme
     * {@link #computeScoreOfConfiguration(List, Set, int)}.
     *
     * @param longReads             long read alignments
     * @param header                header for the long reads
     * @param scoreDiffTolerance    a tolerance where if two configurations' scores differ less than or equal to this amount, they are considered equally good
     * @param toolLogger            logger for, most likely, debugging uses
     *
     * @return              contigs with alignments filtered and custom formatted as {@link AlignmentInterval}
     */
    static JavaRDD<AlignedContig> filterByScore(final JavaRDD<GATKRead> longReads,
                                                final SAMFileHeader header,
                                                final String nonCanonicalContigNamesFile,
                                                final Double scoreDiffTolerance,
                                                final Logger toolLogger) {

        longReads.cache();
        toolLogger.info( "Processing " + longReads.count() + " raw alignments from " +
                         longReads.map(GATKRead::getName).distinct().count() + " contigs.");

        final JavaRDD<AlignedContig> parsedContigAlignments =
                new SvDiscoverFromLocalAssemblyContigAlignmentsSpark.SAMFormattedContigAlignmentParser(longReads, header, false)
                        .getAlignedContigs()
                        .filter(FilterLongReadAlignmentsSAMSpark::contigFilter).cache();
        longReads.unpersist();
        toolLogger.info( "Primitive filtering based purely on MQ left " + parsedContigAlignments.count() + " contigs.");

        return filterAndSplitGappedAI(parsedContigAlignments, nonCanonicalContigNamesFile, header.getSequenceDictionary(), scoreDiffTolerance);
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
                                                         final SAMSequenceDictionary dictionary,
                                                         final Double scoreDiffTolerance) {

        final Set<String> canonicalChromosomes = getCanonicalChromosomes(nonCanonicalContigNamesFile, dictionary);

        return parsedContigAlignments
                .mapToPair(alignedContig -> new Tuple2<>(alignedContig.contigName,
                        new Tuple2<>(alignedContig.contigSequence, pickBestConfigurations(alignedContig, canonicalChromosomes, scoreDiffTolerance))))
                .flatMap(FilterLongReadAlignmentsSAMSpark::reConstructContigFromPickedConfiguration);
    }

    /**
     * Pick the best configurations based on a heuristic scoring scheme implemented in
     * {@link #computeScoreOfConfiguration(List, Set, int)}.
     * @return a 2-D list, where in the case when multiple configurations are equally top-scored, all such configurations are picked up
     */
    @VisibleForTesting
    public static List<List<AlignmentInterval>> pickBestConfigurations(final AlignedContig alignedContig,
                                                                final Set<String> canonicalChromosomes,
                                                                final Double scoreDiffTolerance) {
        return pickBestConfigurations(alignedContig.contigName, alignedContig.alignmentIntervals, canonicalChromosomes, scoreDiffTolerance);
    }

    public static List<List<AlignmentInterval>> pickBestConfigurations(final String contigName,
                                                                       final List<AlignmentInterval> intervals,
        final Set<String> canonicalChromosomes, final Double scoreDiffTolerance) {


            // group 1: get max aligner score of mappings to canonical chromosomes and speed up in case of too many mappings
        final int maxCanonicalChrAlignerScore = intervals.stream()
                .filter(alignmentInterval -> canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig()))
                .mapToInt(ai -> ai.alnScore).max().orElse(0); // possible that no mapping to canonical chromosomes

        // speed up if number of alignments is too high (>10)
        // if mapped to canonical chromosomes, MQ must be >10; otherwise, must have AS higher than max canonical aligner score
        final List<AlignmentInterval> alignmentIntervals;
        if (intervals.size() > 10) {
            alignmentIntervals = intervals.stream()
                    .filter(alignmentInterval -> (!canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig())
                                                                           && alignmentInterval.alnScore > maxCanonicalChrAlignerScore)
                                                 || alignmentInterval.mapQual>10)
                    .collect(Collectors.toList());
        } else {
            alignmentIntervals = intervals;
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
        final Comparator<AlignedContig> numFirst
                = (AlignedContig x, AlignedContig y) -> Integer.compare(x.alignmentIntervals.size(), y.alignmentIntervals.size());
        final Comparator<AlignedContig> mismatchSecond
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

            try (final Stream<String> nonCanonical = Files.lines(IOUtils.getPath((nonCanonicalContigNamesFile)))) {
                return new HashSet<>( Sets.difference(dictionary.getSequences().stream().map(SAMSequenceRecord::getSequenceName).collect(Collectors.toSet()),
                        nonCanonical.collect(Collectors.toSet())) );
            } catch ( final IOException ioe ) {
                throw new UserException("Can't read nonCanonicalContigNamesFile file "+nonCanonicalContigNamesFile, ioe);
            }
        } else {
            final List<String> first22ChromosomesNum = IntStream.range(0, 23).mapToObj(String::valueOf).collect(Collectors.toList());
            final Set<String> canonicalChromosomeNames = first22ChromosomesNum.stream().map(name -> "chr" + name).collect(Collectors.toSet());
            canonicalChromosomeNames.addAll(first22ChromosomesNum);
            canonicalChromosomeNames.addAll(Arrays.asList("chrX", "chrY", "chrM", "X", "Y", "MT"));
            return new HashSet<>( canonicalChromosomeNames );
        }
    }
}
