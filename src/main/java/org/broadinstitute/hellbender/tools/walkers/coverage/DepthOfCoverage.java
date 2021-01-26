package org.broadinstitute.hellbender.tools.walkers.coverage;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.refseq.RefSeqFeature;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.IOException;
import java.util.*;

/**
 * Assess sequence coverage by a wide array of metrics, partitioned by sample, read group, or library
 *
 * <p>
 * This tool processes a set of bam files to determine coverage at different levels of partitioning and
 * aggregation. Coverage can be analyzed per locus, per interval, per gene, or in total; can be partitioned by
 * sample, by read group, by technology, by center, or by library; and can be summarized by mean, median, quartiles,
 * and/or percentage of bases covered to or beyond a threshold. Additionally, reads and bases can be filtered by
 * mapping or base quality score.
 * </p>
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>One or more bam files (with proper headers) to be analyzed for coverage statistics</li>
 *     <li>(Optional) A REFSEQ file to aggregate coverage to the gene level (for information about creating the REFSEQ file, please consult the online documentation)</li>
 * </ul>

 * <h3>Output</h3>
 * <p>
 * Tables pertaining to different coverage summaries. Suffix on the table files declares the contents:
 * </p>
 * <ul>
 *     <li>no suffix: per locus coverage</li>
 *     <li>_summary: total, mean, median, quartiles, and threshold proportions, aggregated over all bases</li>
 *     <li>_statistics: coverage histograms (# locus with X coverage), aggregated over all bases</li>
 *     <li>_interval_summary: total, mean, median, quartiles, and threshold proportions, aggregated per interval</li>
 *     <li>_interval_statistics: 2x2 table of # of intervals covered to >= X depth in >=Y samples</li>
 *     <li>_gene_summary: total, mean, median, quartiles, and threshold proportions, aggregated per gene</li>
 *     <li>_gene_statistics: 2x2 table of # of genes covered to >= X depth in >= Y samples. In its current incarnation it will not include genes not at least partially covered (see --omit-genes-not-entirely-covered-by-traversal for details) </li>
 *     <li>_cumulative_coverage_counts: coverage histograms (# locus with >= X coverage), aggregated over all bases</li>
 *     <li>_cumulative_coverage_proportions: proprotions of loci with >= X coverage, aggregated over all bases</li>
 * </ul>
 *
 * <h3>Notes</h3>
 * <ul>
 *     <li>DepthOfCoverage currently only supports typical nucleotide (and N) bases, IUPAC ambiguity codes or other non-ATCGN bases will cause exceptions</li>
 *     <li>Read filters are applied to the reads before being counted in coverage information. By default Duplicate Marked and non-primary alignments are not counted. This can be disabled with --disable-tool-default-read-filters.</li>
 *     <li>In order to filter reads out by their mapping qualities, the recommended approach is to use the MappingQualityReadFilter with the --minimum-mapping-quality or --maximum-mapping-quality arguments specified</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk \
 *   DepthOfCoverage \
 *   -R reference.fasta \
 *   -O file_name_base \
 *   -I input_bams.list
 *   [-gene-list refSeq.sorted.refseq] \
 *   [-partition-type readgroup] \
 *   [--summary-coverage-threshold 4 --summary-coverage-threshold 6 --summary-coverage-threshold 10] \
 *   [-L my_capture_genes.interval_list]
 * </pre>
 */

@CommandLineProgramProperties(
        summary = "Generate coverage summary information for reads data",
        oneLineSummary = "Generate coverage summary information for reads data",
        programGroup = CoverageAnalysisProgramGroup.class)
@BetaFeature
@DocumentedFeature
public class DepthOfCoverage extends LocusWalkerByInterval {
    private CoverageOutputWriter writer;
    // Map used to store running aggregate counts for intervals that are being recorded by DepthOfCoverage
    private Map<Locatable, DepthOfCoveragePartitionedDataStore> activeCoveragePartitioner = new HashMap<>();
    // Map used to store target tables for each partition which are used in the computation of median coverage scores
    private Map<DoCOutputType.Partition, int[][]> perIntervalStatisticsAggregationByPartitioning = new HashMap<>();
    // Map used to store target tables for each partition which are used in the computation of median coverage scores
    private Map<DoCOutputType.Partition, int[][]> perGeneStatisticsAggregationByPartitioning = new HashMap<>();

    // PartitionDataStore corresponding to every base traversed by the tool
    private DepthOfCoveragePartitionedDataStore coverageTotalsForEntireTraversal;
    // List of all of the samples to be output split by the partition type
    private Map<DoCOutputType.Partition, List<String>> globalIdentifierMap;

    /**
     * Base file name about which to create the coverage information
     */
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Base file location to which to write coverage summary information, must be a path to a file")
    private String baseFileName = null;
    /**
     * Bases with quality scores lower than this threshold will be skipped. This is set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "min-base-quality", doc = "Minimum quality of bases to count towards depth", optional = true, minValue = 0, maxValue = Byte.MAX_VALUE)
    private byte minBaseQuality = 0;
    /**
     * Bases with quality scores higher than this threshold will be skipped. The default value is the largest number that can be represented as a byte.
     */
    @Argument(fullName = "max-base-quality", doc = "Maximum quality of bases to count towards depth", optional = true, minValue = 0, maxValue = Byte.MAX_VALUE)
    private byte maxBaseQuality = Byte.MAX_VALUE;

    @Argument(fullName = "count-type", doc = "How should overlapping reads from the same fragment be handled? NOTE: currently only COUNT_READS is supported", optional = true)
    private CoverageUtils.CountPileupType countType = CoverageUtils.CountPileupType.COUNT_READS;

    /**
     * Instead of reporting depth, the program will report the base pileup at each locus
     */
    @Argument(fullName = "print-base-counts", doc = "Add base counts to per-locus output", optional = true)
    private boolean printBaseCounts = false;

    /**
     * Disabling the tabulation of locus statistics (# loci covered by sample by coverage) should speed up processing.
     */
    @Argument(fullName = "omit-locus-table", doc = "Do not calculate per-sample per-depth counts of loci", optional = true)
    private boolean omitLocusTable = false;

    /**
     * Disabling the tabulation of interval statistics (mean, median, quartiles AND # intervals by sample by coverage) should speed up processing.
     */
    @Argument(fullName = "omit-interval-statistics", doc = "Do not calculate per-interval statistics", optional = true, mutex = "calculate-coverage-over-genes")
    private boolean omitIntervals = false;
    /**
     * Disabling the tabulation of total coverage at every base should speed up processing.
     */
    @Argument(fullName = "omit-depth-output-at-each-base", doc = "Do not output depth of coverage at each base", optional = true)
    private boolean omitDepthOutput = false;
    /**
     * This option simply disables writing separate files for per-sample summary statistics (total, mean, median, quartile coverage per sample). These statistics are still calculated internally, so enabling this option will not improve runtime.
     */
    @Argument(fullName = "omit-per-sample-statistics", doc = "Do not output the summary files per-sample", optional = true)
    private boolean omitSampleSummary = false;

    /**
     * Specify a RefSeq file for use in aggregating coverage statistics over genes.
     * <p>
     * This argument is incompatible with --omit-interval-statistics. A warning will be logged and no output file will be produced for the gene list if these arguments are enabled together.
     */
    @Argument(fullName = "calculate-coverage-over-genes", shortName = "gene-list", doc = "Calculate coverage statistics over this list of genes", optional = true, mutex = "omit-interval-statistics")
    private List<String> refSeqGeneListFiles = new ArrayList<>();
    /**
     * Remove genes from the gene summary output file if all of its exon bases were not completely covered by traversal.
     */
    @Argument(fullName = "omit-genes-not-entirely-covered-by-traversal", doc = "Do not output gene summary if it was not completely covered by traversal intervals", optional = true)
    private boolean omitPartiallyCoveredGenes = false;

    /**
     * Output file format (e.g. csv, table, rtable); defaults to r-readable table.
     */
    @Argument(fullName = "output-format", doc = "The format of the output file", optional = true)
    CoverageOutputWriter.DEPTH_OF_COVERAGE_OUTPUT_FORMAT outputFormat = CoverageOutputWriter.DEPTH_OF_COVERAGE_OUTPUT_FORMAT.CSV;


    // ---------------------------------------------------------------------------
    //
    // Advanced arguments
    //
    // ---------------------------------------------------------------------------

    /**
     * Normally, sites where the reference is N (or another non-canonical base) are skipped. If this option is enabled, these sites will be included in DoC calculations if there is coverage from neighboring reads.
     */
    @Advanced
    @Argument(fullName = "include-ref-n-sites", doc = "Include sites where the reference is N", optional = true)
    private boolean includeRefNBases = false;
    /**
     * Sets the low-coverage cutoff for granular binning. All loci with depth < START are counted in the first bin.
     */
    @Advanced
    @Argument(fullName = "start", doc = "Starting (left endpoint) for granular binning", optional = true, minValue = 0)
    private int start = 1;
    /**
     * Sets the high-coverage cutoff for granular binning. All loci with depth > STOP are counted in the last bin.
     */
    @Advanced
    @Argument(fullName = "stop", doc = "Ending (right endpoint) for granular binning", optional = true, minValue = 1)
    private int stop = 500;
    /**
     * Sets the number of bins for granular binning
     */
    @Advanced
    @Argument(fullName = "nBins", doc = "Number of bins to use for granular binning", optional = true, minValue = 0, minRecommendedValue = 1)
    private int nBins = 499;

    /**
     * By default, coverage is partitioned by sample, but it can be any combination of sample, readgroup and/or library.
     */
    @Argument(fullName = "partition-type", shortName = "pt", doc = "Partition type for depth of coverage", optional = true)
    private EnumSet<DoCOutputType.Partition> partitionTypes = EnumSet.of(DoCOutputType.Partition.sample);

    /**
     * Consider a spanning deletion as contributing to coverage. Also enables deletion counts in per-base output.
     */
    @Advanced
    @Argument(fullName = "include-deletions", doc = "Include information on deletions alongside other bases in output table counts", optional = true)
    private boolean includeDeletions = false;

    @Advanced
    @Argument(fullName = "ignore-deletion-sites", doc = "Ignore sites consisting only of deletions", optional = true)
    boolean ignoreDeletionSites = false;

    /**
     * For summary file outputs, report the percentage of bases covered to an amount equal to or greater than this number  (e.g. % bases >= CT for each sample). Defaults to 15; can take multiple arguments.
     */
    @Advanced
    @Argument(fullName = "summary-coverage-threshold", doc = "Coverage threshold (in percent) for summarizing statistics", optional = true)
    private List<Integer> coverageThresholds = new ArrayList<>(Collections.singleton(15));

    // We explicitly handle read N bases in the code
    @Override
    public boolean includeNs() {
        return true;
    }

    // We want to make sure to still generate coverage information over uncovered bases.
    @Override
    public boolean emitEmptyLoci() {
        return true;
    }

    @Override
    public boolean includeDeletions() {
        return includeDeletions && ! ignoreDeletionSites;
    }

    // For now this tool requires a reference: s
    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> defaultFilters = new ArrayList<>(2);
        defaultFilters.add(new WellformedReadFilter());
        defaultFilters.add(new ReadFilterLibrary.NotDuplicateReadFilter());
        defaultFilters.add(ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT);
        defaultFilters.add(new ReadFilterLibrary.MappedReadFilter());
        return defaultFilters;
    }


    @Override
    public void onTraversalStart() {

        try {
            writer = new CoverageOutputWriter(outputFormat,
                    partitionTypes,
                    baseFileName,
                    !refSeqGeneListFiles.isEmpty(),
                    printBaseCounts,
                    omitDepthOutput,
                    omitIntervals,
                    omitSampleSummary,
                    omitLocusTable,
                    coverageThresholds);
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Couldn't create output file, encountered exception: " + e.getMessage(), e);
        }

        globalIdentifierMap = makeGlobalIdentifierMap(partitionTypes);
        writer.writeCoverageOutputHeaders(globalIdentifierMap);

        ReadUtils.getSamplesFromHeader(getHeaderForReads());

        coverageTotalsForEntireTraversal = createCoveragePartitioner();
    }

    /**
     *
     * @param alignmentContext current alignment context
     * @param referenceContext Reference bases spanning the current locus. Will be an empty, but non-null, context object
     *                         if there is no backing source of reference data (in which case all queries on it will return
     *                         an empty array/iterator). Can request extra bases of context around the current locus
     *                         by invoking {@link ReferenceContext#setWindow} on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext Features spanning the current locus. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     * @param activeIntervals
     */
    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext, Set<Locatable> activeIntervals) {
        // TODO evaluate consequences of supporting nonexistant references
        if (includeRefNBases || (hasReference() && BaseUtils.isRegularBase(referenceContext.getBase()))) {
            final Map<DoCOutputType.Partition, Map<String, int[]>> countsByPartition = CoverageUtils.getBaseCountsByPartition(alignmentContext, minBaseQuality, maxBaseQuality, countType, partitionTypes, getHeaderForReads());

            if (!omitDepthOutput) {
                writer.writePerLocusDepthSummary(referenceContext.getInterval(), countsByPartition, globalIdentifierMap, includeDeletions);
            }

            // Update the traversing partitioners with this locus data:
            coverageTotalsForEntireTraversal.addLocusData(countsByPartition);

            // Update all of the active intervals that we are tracking seperately with the generated counts
            for (Locatable loc : activeCoveragePartitioner.keySet()) {
                // For genes, we don't want to update the interval for non-exon bases
                if (loc.contains(alignmentContext)) {
                    activeCoveragePartitioner.get(loc).addLocusData(countsByPartition);
                }
            }
        }
    }

    /**
     * Ensure that we are tracking the newly provided interval with per-locus data
     *
     * @param activeInterval Locatable provided to the walker to be initialized
     */
    public void onIntervalStart(Locatable activeInterval) {
        DepthOfCoveragePartitionedDataStore partitioner = createCoveragePartitioner();
        //NOTE: we don't populate perIntervalStatisticsAggregationByPartitioning here because that gets populated on the fly later
        activeCoveragePartitioner.put(activeInterval, partitioner);
    }

    /**
     * When this method is called, we expect at some point in the past onIntervalStart() was called on the same Locatable
     * object and that we have been continuously tracking and adding per-base information for each base in the pileup
     * that overlaps those objects. Given the interval being closed out, we perform one of the following tasks:
     *
     *  - Write out per-interval summary information to each of the partition interval summary files and update the target
     *    summary table corresponding to each track.
     *
     *  - If it corresponded to a gene interval, write out to the gene summary tables the line corresponding the coverage
     *    over targets contained within the genes exome.
     *
     * @param activeInterval
     */
    public void onIntervalEnd(final Locatable activeInterval) {
        final DepthOfCoveragePartitionedDataStore partitionerToRemove = activeCoveragePartitioner.remove(activeInterval);
        if (activeInterval instanceof SimpleInterval) {
            // For each partition type we are managing, make sure we update the active statistics table
            if (!omitIntervals) {
                for (DoCOutputType.Partition p : partitionTypes) {
                    // Write the per-interval depth information as necessary
                    final DepthOfCoverageStats coverageByAggregationPartitionType = partitionerToRemove.getCoverageByAggregationType(p);
                    writer.writePerIntervalDepthInformation(p, (SimpleInterval) activeInterval, coverageByAggregationPartitionType, globalIdentifierMap.get(p));

                    // Create a new table if necessary
                    if (!perIntervalStatisticsAggregationByPartitioning.containsKey(p)) {
                        perIntervalStatisticsAggregationByPartitioning.put(p, new int[coverageByAggregationPartitionType.getHistograms().size()][coverageByAggregationPartitionType.getEndpoints().length + 1]);
                    }
                    // Update the target table to reflect the updated coverage information for this target
                    CoverageUtils.updateTargetTable(perIntervalStatisticsAggregationByPartitioning.get(p), coverageByAggregationPartitionType);
                }
            }

        } else if (activeInterval instanceof RefSeqFeature) {
            DepthOfCoverageStats coverageBySample = partitionerToRemove.getCoverageByAggregationType(DoCOutputType.Partition.sample);

            if ( ! omitPartiallyCoveredGenes || ((RefSeqFeature)activeInterval).getTotalExonLength() <= coverageBySample.getNumberOfLociCovered()) {
                writer.writePerGeneDepthInformation((RefSeqFeature) activeInterval, coverageBySample, globalIdentifierMap.get(DoCOutputType.Partition.sample));

                final DepthOfCoverageStats coverageByAggregationPartitionType = partitionerToRemove.getCoverageByAggregationType(DoCOutputType.Partition.sample);

                // Create a new table if necessary
                if (!perGeneStatisticsAggregationByPartitioning.containsKey(DoCOutputType.Partition.sample)) {
                    perGeneStatisticsAggregationByPartitioning.put(DoCOutputType.Partition.sample, new int[coverageByAggregationPartitionType.getHistograms().size()][coverageByAggregationPartitionType.getEndpoints().length + 1]);
                }

                // Update the target table to reflect the updated coverage information for this target
                CoverageUtils.updateTargetTable(perGeneStatisticsAggregationByPartitioning.get(DoCOutputType.Partition.sample), coverageByAggregationPartitionType);

            }
        } else {
            throw new GATKException("Unrecognized Locatable object supplied for traversal, only RefSeqFeature and SimpleInterval are supported: "+activeInterval.toString());
        }
    }

    /**
     * When we finish traversal we should have the following held in state:
     *  1. {@link #perIntervalStatisticsAggregationByPartitioning} that contains a map from partition to a doubly indexed integer array
     *     that corresponds to a count of the number of data sources with at least the given depth. Importantly this is
     *     expected to have had CoverageUtils.updateTargetTable() called on it once per interval (i.e. once per locus).
     *
     *  2. {@link #coverageTotalsForEntireTraversal} Which we expect to have had update() called on it exactly once for
     *     every locus traversed by the locus walker. Consequently we expect partitioner to contain partitioned counts
     *     corresponding to the entire window covered by the tool.
     *
     */
    @Override
    public Object onTraversalSuccess() {
        // Write out accumulated interval summary statistics
        for (DoCOutputType.Partition partition : perIntervalStatisticsAggregationByPartitioning.keySet()) {
            writer.writeOutputIntervalStatistics(partition, perIntervalStatisticsAggregationByPartitioning.get(partition), CoverageUtils.calculateCoverageHistogramBinEndpoints(start,stop,nBins));
        }

        // Write out accumulated gene summary statistics
        if (!refSeqGeneListFiles.isEmpty()) {
            for (DoCOutputType.Partition partition : perGeneStatisticsAggregationByPartitioning.keySet()) {
                writer.writeOutputGeneStatistics(perGeneStatisticsAggregationByPartitioning.get(partition), CoverageUtils.calculateCoverageHistogramBinEndpoints(start, stop, nBins));
            }
        }

        if (!omitSampleSummary) {
            logger.info("Outputting summary info");
            for (DoCOutputType.Partition type : partitionTypes) {
                writer.writeCumulativeOutputSummaryFiles(coverageTotalsForEntireTraversal, type, globalIdentifierMap.get(type));
            }
        }

        if (!omitLocusTable) {
            logger.info("Outputting locus summary");
            for (DoCOutputType.Partition type : partitionTypes) {
                writer.writePerLocusCumulativeCoverageMetrics(coverageTotalsForEntireTraversal, type, globalIdentifierMap.get(type));
            }
        }

        writer.close();
        return "success";
    }

    // Initialize a coveragePartitioner object for storing interval data.
    private DepthOfCoveragePartitionedDataStore createCoveragePartitioner() {
        return new DepthOfCoveragePartitionedDataStore(partitionTypes, start, stop, nBins, includeDeletions, omitLocusTable, globalIdentifierMap);
    }

    /**
     * This makes a global sample map by identifier. That is, for every requested partition all of the samples corresponding
     * to that partition are listed and sorted lexicographically.
     *
     * @param types Partition types to generate samples for
     * @return A map of {@link DoCOutputType.Partition} to list of sorted samples
     */
    private Map<DoCOutputType.Partition, List<String>> makeGlobalIdentifierMap(Collection<DoCOutputType.Partition> types) {
        Map<DoCOutputType.Partition, List<String>> partitions = new LinkedHashMap<>();
        for (DoCOutputType.Partition t : types) {
            List<String> samplesForType = new ArrayList<>(getSamplesByPartitionFromReadHeader(t));
            samplesForType.sort(String::compareTo);
            partitions.put(t, Collections.unmodifiableList(samplesForType));
        }
        return Collections.unmodifiableMap(partitions);
    }

    // Return all of the samples denoted by a given partition based on the header.
    private HashSet<String> getSamplesByPartitionFromReadHeader(DoCOutputType.Partition type) {
        final HashSet<String> partition = new HashSet<String>();
        final SAMFileHeader header = getHeaderForReads();
        if (type == DoCOutputType.Partition.sample) {
            partition.addAll(ReadUtils.getSamplesFromHeader(header));
        } else {
            for (SAMReadGroupRecord rg : header.getReadGroups()) {
                partition.add(CoverageUtils.getTypeID(rg, type));
            }
        }
        return partition;
    }


    /**
     * We want to keep track of per-interval information for both the user specified intervals and for user specified genes
     *
     * @return A combined list of the unmerged user specified intervals and any specified refSeqFeatures specified
     */
    @Override
    public List<Locatable> getIntervalObjectsToQueryOver() {
        List<? extends Locatable> userProvidedIntervals = intervalArgumentCollection.getIntervalsWithoutMerging(getBestAvailableSequenceDictionary());

        final List<Locatable> refSeqInputs = new ArrayList<>(userProvidedIntervals);
        for (String input : refSeqGeneListFiles) {
            FeatureDataSource<RefSeqFeature> refSeqReader = new FeatureDataSource<>(input);
            for (final RefSeqFeature geneRecord : refSeqReader) {
                refSeqInputs.add(geneRecord);
            }
            refSeqReader.close();
        }
        return refSeqInputs;
    }
}


