package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.ChimericAlignment;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.*;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoverFromLocalAssemblyContigAlignmentsSpark.SAMFormattedContigAlignmentParser;


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
@ExperimentalFeature
@CommandLineProgramProperties(
        oneLineSummary = "Examines alignments of chimeric contigs, attempting to produce an optimal tiling",
        summary =
        "This tool is used in development and should not be of interest to most researchers.  It is a prototype" +
        " of one aspect of structural variant calling from chimeric contigs produced by local assemblies." +
        " It takes a SAM/BAM/CRAM containing the alignments of assembled contigs and filters them with the" +
        " aim of providing \"optimal coverage\" of the contig, based on an heuristic scoring scheme.",
        programGroup = DiagnosticsAndQCProgramGroup.class)
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
                    () -> AssemblyContigAlignmentsConfigPicker
                            .createOptimalCoverageAlignmentSetsForContigs(reads, header, nonCanonicalContigNamesFile, configScoreDiffTolerance, localLogger)
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
                new SAMFormattedContigAlignmentParser(longReads, header, true)
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


}
