package org.broadinstitute.hellbender.tools.walkers.coverage;

import com.google.common.collect.Lists;
import com.opencsv.CSVWriter;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.refseq.RefSeqFeature;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.*;
import java.math.RoundingMode;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.*;


/**
 * This is a class for managing the output formatting/files for DepthOfCoverage.
 *
 * Output for {@link DepthOfCoverage} can be organized as a combination Partition, Aggregation, and OutputType, with this writer
 * storing an internal list of streams to which to write output organized by {@link DoCOutputType} objects. Generally
 * speaking there are three patterns that this writer is responsible for managing:
 *
 *  1. writePerLocusDepthSummary() - This should be called over every locus and is responsible for producing the locus
 *                                   coverage output table summarizing the coverage (and possibly base counts)
 *                                   for every partition type x relevant samples. This takes as input the output from
 *                                   {@link CoverageUtils#getBaseCountsByPartition}
 *
 *  2a. writePerIntervalDepthInformation() - This should be called once per traversal interval once it is finished and takes
 *                                        as input a {@link DepthOfCoveragePartitionedDataStore} object corresponding to the coverage
 *                                        information over the whole interval. This is used to write both "_interval_summary"
 *                                        and "_interval_statistics" files for each partition.
 *
 *  2b. writePerGeneDepthInformation() - Similarly to 2a, this should be called once per gene in the interval traversal and
 *                                       it takes the same gene-interval coverage summary as 2a.
 *                                       NOTE: This may not actually be called in function for every gene if the provided traversal
 *                                             intervals for the tool do not actually cover any bases for the gene in quesion. If
 *                                             this is the case then the gene will be passed over. //TODO in the future this should be changed to emit empty coverage counts for untraversed genes.
 *
 *  3. writeTraversalCumulativeCoverage() - This method takes a {@link DepthOfCoveragePartitionedDataStore} object that should correspond
 *                                         to the partitioned counts for every base traversed by DepthOfCoverage aggregated.
 *                                         This method is responsible for outputting the "_cumulative_coverage_counts",
 *                                         "_cumulative_coverage_proportions", "_statistics", and "_summary" files.
 *
 */
public class CoverageOutputWriter implements Closeable {
    static final List<String> CUMULATIVE_SUMMARY_OUTPUT_LINES = Collections.unmodifiableList(Arrays.asList("sample_id", "total", "mean", "granular_third_quartile", "granular_median", "granular_first_quartile"));
    private static final char TSV_SEPARATOR = '\t';
    private static final char CSV_SEPARATOR = ',';

    private final EnumSet<DoCOutputType.Partition> partitions;
    private final boolean includeGeneOutput;
    private final boolean omitIntervals;
    private List<String> summaryHeaderSampleSuffixes;
    private boolean printBaseCounts;
    private Map<DoCOutputType, SimpleCSVWriterWrapperWithHeader> outputs;
    private DEPTH_OF_COVERAGE_OUTPUT_FORMAT outputFormat;
    private boolean omitDepthOutput;
    private List<Integer> coverageThresholds;
    private final DecimalFormat DOUBLE_FORMAT_2PLACES = new DecimalFormat("0.00");
    private final DecimalFormat DOUBLE_FORMAT_1PLACE = new DecimalFormat("0.0");

    public enum DEPTH_OF_COVERAGE_OUTPUT_FORMAT {
        TABLE,
        CSV
    }


    /**
     * Creates a CoverageOutputWriter for managing all of the output files for DepthOfCoverage.
     * <p>
     * This object holds on to output writers for each of the tables and ensures, among other things, that the data are
     * written into the correct table corresponding to each data output call. This class is also responsible for generating
     * the data needed to format some output lines for DepthOfCoverage. This means that for most datatypes, DepthOfCoverage
     * need only provided with the proper {@link DepthOfCoverageStats} object and the this class will handle extracting
     * the necessary summary data therin. Furthermore this class understands how to write CSV/TSV files, which all DoC output
     * files approximately follow.
     *
     * @param outputFormat               type of table (CSV/TSV) to output results as
     * @param partitionsToCover          partitioning for the data (eg sample, library, etc...) that will be computed in paralell
     * @param outputBaseName             base file path for the output tables (will be prepended to the output suffix and .tsv or .csv)
     * @param includeGeneOutputPerSample whether to produce an output sink for gene partition data
     * @param printBaseCounts            whether to summarize per-nucleotide counts at each locus
     * @param omitDepthOutput            if true will not generate locus output files
     * @param omitIntervals              if true then will generate "_interval_summary" output files
     * @param omitSampleSummary          if true then will not generate "_statistics" and "_summary" output sinks
     * @param omitLocusTable             if true will not generate "_cumulative_coverage_counts" or "_cumulative_coverage_proportions"
     * @param coverageThresholds         Threshold coverage level to use for reporting.
     * @throws IOException
     */
    public CoverageOutputWriter(final DEPTH_OF_COVERAGE_OUTPUT_FORMAT outputFormat,
                                final EnumSet<DoCOutputType.Partition> partitionsToCover,
                                final String outputBaseName,
                                final boolean includeGeneOutputPerSample,
                                final boolean printBaseCounts,
                                final boolean omitDepthOutput,
                                final boolean omitIntervals,
                                final boolean omitSampleSummary,
                                final boolean omitLocusTable,
                                final List<Integer> coverageThresholds) throws IOException {
        // NOTE: we have to set the rounding mode here in order to preserve rounding behavior to match test results based on GATK3 which uses HALF_UP rounding
        DOUBLE_FORMAT_2PLACES.setRoundingMode(RoundingMode.HALF_UP);
        DOUBLE_FORMAT_1PLACE.setRoundingMode(RoundingMode.HALF_UP);
        this.outputFormat = outputFormat;
        this.partitions = partitionsToCover;
        this.includeGeneOutput = includeGeneOutputPerSample;
        this.omitIntervals = omitIntervals;
        this.printBaseCounts = printBaseCounts;
        this.omitDepthOutput = omitDepthOutput;
        this.coverageThresholds = coverageThresholds;
        this.summaryHeaderSampleSuffixes = getSampleSuffixes(coverageThresholds);

        char separator;
        if (outputFormat == DEPTH_OF_COVERAGE_OUTPUT_FORMAT.CSV) {
            separator = CSV_SEPARATOR;
        } else {
            separator = TSV_SEPARATOR;
        }

        outputs = new HashMap<>();
        if (!omitDepthOutput) {
            // Create a depth summary output sink
            DoCOutputType depthSummaryByLocus = new DoCOutputType(null, DoCOutputType.Aggregation.locus, DoCOutputType.FileType.summary);
            outputs.put(depthSummaryByLocus, getOutputStream(outputBaseName, depthSummaryByLocus, separator));
        }

        if (!omitIntervals) {
            for (DoCOutputType.Partition partition : partitionsToCover) {
                // Create an interval summary output sink
                DoCOutputType intervalSummaryBypartition = new DoCOutputType(partition, DoCOutputType.Aggregation.interval, DoCOutputType.FileType.summary);
                outputs.put(intervalSummaryBypartition, getOutputStream(outputBaseName, intervalSummaryBypartition, separator));

                // Create an interval statistics output sink
                DoCOutputType intervalStatisticsByPartition = new DoCOutputType(partition, DoCOutputType.Aggregation.interval, DoCOutputType.FileType.statistics);
                outputs.put(intervalStatisticsByPartition, getOutputStream(outputBaseName, intervalStatisticsByPartition, separator));
            }
        }

        // The gene output summary is computed on a per-sample bases and is related to the per-sample output tables. Thus we require sample partitioning be present.
        if (canWriteGeneOutput()) {
            // Create a special output sink for gene data if provided
            DoCOutputType geneSummaryOut = new DoCOutputType(DoCOutputType.Partition.sample, DoCOutputType.Aggregation.gene, DoCOutputType.FileType.summary);
            outputs.put(geneSummaryOut, getOutputStream(outputBaseName, geneSummaryOut, separator));

            // Create an gene statistics output sink
            DoCOutputType geneStatisticsByPartition = new DoCOutputType(DoCOutputType.Partition.sample, DoCOutputType.Aggregation.gene, DoCOutputType.FileType.statistics);
            outputs.put(geneStatisticsByPartition, getOutputStream(outputBaseName, geneStatisticsByPartition, separator));
        }

        if (!omitSampleSummary) {
            for (DoCOutputType.Partition partition : partitionsToCover) {
                // Create a cumulative coverage summary output sink
                DoCOutputType cumulativeSummaryOut = new DoCOutputType(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.summary);
                outputs.put(cumulativeSummaryOut, getOutputStream(outputBaseName, cumulativeSummaryOut, separator));

                // Create a cumulative coverage statistics output sink
                DoCOutputType cumulativeStatisticsOut = new DoCOutputType(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.statistics);
                outputs.put(cumulativeStatisticsOut, getOutputStream(outputBaseName, cumulativeStatisticsOut, separator));
            }
        }

        if (!omitLocusTable) {
            for (DoCOutputType.Partition partition : partitionsToCover) {
                // Create a cumulative coverage counts output sink
                DoCOutputType cumulativeCoverageCountsOut = new DoCOutputType(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.coverage_counts);
                outputs.put(cumulativeCoverageCountsOut, getOutputStream(outputBaseName, cumulativeCoverageCountsOut, separator));

                // Create a cumulative coverage portions output sink
                DoCOutputType cumulativeCoverageProportionsOut = new DoCOutputType(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.coverage_proportions);
                outputs.put(cumulativeCoverageProportionsOut, getOutputStream(outputBaseName, cumulativeCoverageProportionsOut, separator));
            }
        }
    }

    // A helper method that determines if we should be handling gene output information.
    private boolean canWriteGeneOutput() {
        return includeGeneOutput && partitions.contains(DoCOutputType.Partition.sample);
    }

    // Helper method to generate an output stream given a DoCOutputType Object and the base name
    private static SimpleCSVWriterWrapperWithHeader getOutputStream(String outputBaseName, DoCOutputType depthSummaryByLocus, char separator) throws IOException {
        return new SimpleCSVWriterWrapperWithHeader(Files.newBufferedWriter((IOUtils.getPath(depthSummaryByLocus.getFilePath(outputBaseName)))), separator);
    }

    /**
     * Initialize and write the output headers for the output streams that need to be continuously updated during traversal.
     *
     * @param sortedSamplesByPartition global map of partitionType to sorted list of samples associated with that partition
     */
    public void writeCoverageOutputHeaders(final Map<DoCOutputType.Partition, List<String>> sortedSamplesByPartition) {
        if (!omitDepthOutput) { // print header
            writePerLocusDepthOutputSummaryHeader(sortedSamplesByPartition);
        }

        // write
        if (canWriteGeneOutput()) {
            final SimpleCSVWriterWrapperWithHeader geneSummaryOut = getCorrectOutputWriter(DoCOutputType.Partition.sample, DoCOutputType.Aggregation.gene, DoCOutputType.FileType.summary);
            geneSummaryOut.addHeaderLine(getIntervalSummaryHeader("Gene", sortedSamplesByPartition.get(DoCOutputType.Partition.sample)));
        }

        if (!omitIntervals) {
            for (DoCOutputType.Partition partition : partitions) {
                final SimpleCSVWriterWrapperWithHeader intervalSummaryOut = getCorrectOutputWriter(partition, DoCOutputType.Aggregation.interval, DoCOutputType.FileType.summary);
                intervalSummaryOut.addHeaderLine(getIntervalSummaryHeader("Target", sortedSamplesByPartition.get(partition)));
            }
        }
    }

    /**
     * Add to the per-locus depth summary output writer a header line describing all the samples in the file optionally
     * with extra columns for each sample to store individual pileup basecounts as well.
     *
     * @param identifiersByType global map of partitions to a sorted list of their samples
     */
    private void writePerLocusDepthOutputSummaryHeader(final Map<DoCOutputType.Partition, List<String>> identifiersByType) {
        SimpleCSVWriterWrapperWithHeader out = getCorrectOutputWriter(null, DoCOutputType.Aggregation.locus, DoCOutputType.FileType.summary);
        List<String> columns = Lists.newArrayList("Locus", "Total_Depth");
        for (DoCOutputType.Partition type : partitions) {
            columns.add("Average_Depth_" + type.toString());
        }

        for (DoCOutputType.Partition type : partitions) {
            for (String s : identifiersByType.get(type)) {
                columns.add("Depth_for_" + s);
                if (printBaseCounts) {
                    columns.add(s + "_base_counts");
                }
            }
        }
        out.addHeaderLine(columns);
    }

    private SimpleCSVWriterWrapperWithHeader getCorrectOutputWriter(DoCOutputType.Partition partition, DoCOutputType.Aggregation aggregation, DoCOutputType.FileType fileType) {
        DoCOutputType outputType = new DoCOutputType(partition, aggregation, fileType);
        if (!outputs.containsKey(outputType)) {
            throw new GATKException(String.format("Unable to find appropriate stream for partition = %s, aggregation = %s, file type = %s", partition, aggregation, fileType));
        }
        return outputs.get(outputType);
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Outward facing writer methods
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Writes a summary of the given per-locus data out to the main output track. Output data for any given sample may
     * optionally contain a summary of which bases were present in the pileup at that site which will consist of a 4
     * element list in the format 'A:0 C:0 T:0 G:0 N:0'.
     * <p>
     * This writer is responsible for determining the total depth at the site by polling one of the output tracks and using
     * that value to calculate the mean per-sample coverage for each partitioning. These depth summary columns as well as
     * the locus site are present in the first lines of the locusSummary output.
     *
     * @param locus                The site corresponding to the data that will be written out
     * @param countsBySampleByType A map of partition to sample and finally to 4 element array of bases. This output is
     *                             expected to be produced by {@link CoverageUtils#getBaseCountsByPartition(AlignmentContext, byte, byte, CoverageUtils.CountPileupType, Collection, SAMFileHeader)}
     * @param identifiersByType    A global map of sorted samples in each partition, to be used for ordering.
     *                             NOTE: this map should remain unchanged between every call of this method or correct output is not guarinteed.
     * @param includeDeletions     Whether or not to include deletions in the summary line
     */
    public void writePerLocusDepthSummary(final SimpleInterval locus, final Map<DoCOutputType.Partition, Map<String, int[]>> countsBySampleByType,
                                          final Map<DoCOutputType.Partition, List<String>> identifiersByType, final boolean includeDeletions) {

        SimpleCSVWriterWrapperWithHeader lineWriter = getCorrectOutputWriter(null, DoCOutputType.Aggregation.locus, DoCOutputType.FileType.summary);
        SimpleCSVWriterWrapperWithHeader.SimpleCSVWriterLineBuilder lineBuilder = lineWriter.getNewLineBuilder();

        // get the depths per sample and buildAndWriteLine up the output string while tabulating total and average coverage
        int tDepth = 0;
        boolean depthCounted = false;
        for (DoCOutputType.Partition type : partitions) {
            Map<String, int[]> countsByID = countsBySampleByType.get(type);
            for (String sample : identifiersByType.get(type)) {
                long dp = (countsByID != null && countsByID.keySet().contains(sample)) ? MathUtils.sum(countsByID.get(sample)) : 0;
                lineBuilder.setColumn("Depth_for_" + sample, Long.toString(dp));
                if (printBaseCounts) {
                    lineBuilder.setColumn(sample + "_base_counts", getBaseCountsString(countsByID != null ? countsByID.get(sample) : null, includeDeletions));
                }
                if (!depthCounted) {
                    tDepth += dp;
                }
            }
            depthCounted = true; // only sum the total depth once
        }

        // Add the total depth and summary information to the line
        lineBuilder.setColumn(0, locus.getContig() + ":" + locus.getStart()).setColumn(1, Integer.toString(tDepth));
        for (DoCOutputType.Partition type : partitions) { //Note that this is a deterministic traversal since the underlying set is an EnumSet
            lineBuilder.setColumn("Average_Depth_" + type.toString(), DOUBLE_FORMAT_2PLACES.format( (double) tDepth / identifiersByType.get(type).size()));
        }
        lineBuilder.buildAndWriteLine();
    }

    /**
     * Method that should be called once per-partition at the end of each coverage interval. This method is responsible
     * for extending the per-interval depth summary information for each sample.
     *
     * @param partition     Partition corresponding to the data to be written
     * @param interval      interval spanned by the input data
     * @param intervalStats statistics for coverage overlapped by the interval
     * @param sortedSamples sorted list of samples for this partition
     */
    public void writePerIntervalDepthInformation(final DoCOutputType.Partition partition, final SimpleInterval interval, final DepthOfCoverageStats intervalStats, final List<String> sortedSamples) {
        SimpleCSVWriterWrapperWithHeader output = getCorrectOutputWriter(partition, DoCOutputType.Aggregation.interval, DoCOutputType.FileType.summary);
        printIntervalSummaryLine(output, interval.toString(), intervalStats, sortedSamples);
    }

    /**
     * Write summary information for a gene. This method should be called for each gene in the coverage input that has been passed in coverage.
     *
     * @param gene          gene corresponding to the coverage statistics
     * @param intervalStats statistics for coverage overlapped by the gene (may or may not be dropped by exons)
     * @param sortedSamples sorted list of samples for this partition
     */
    public void writePerGeneDepthInformation(final RefSeqFeature gene, final DepthOfCoverageStats intervalStats, final List<String> sortedSamples) {
        SimpleCSVWriterWrapperWithHeader output = getCorrectOutputWriter(DoCOutputType.Partition.sample, DoCOutputType.Aggregation.gene, DoCOutputType.FileType.summary);
        printIntervalSummaryLine(output, gene.getGeneName(), intervalStats, sortedSamples);
    }

    /**
     * Write per-locus cumulative coverage metrics. Note that this should be invoked on a coverage partitioner that has
     * been updated for every locus across the traversal intervals.
     *
     * @param coverageProfilesForEntireTraversal DepthOfCoveragePartitionedDataStore object corresponding to the entire tool traversal.
     * @param partition                          Partition corresponding to the data to be written
     * @param sortedSampleLists                  sorted list of samples for this partition
     */
    public void writePerLocusCumulativeCoverageMetrics(final DepthOfCoveragePartitionedDataStore coverageProfilesForEntireTraversal, final DoCOutputType.Partition partition,
                                                       final List<String> sortedSampleLists) {
        outputPerLocusCumulativeSummaryAndStatistics(getCorrectOutputWriter(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.coverage_counts),
                getCorrectOutputWriter(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.coverage_proportions),
                coverageProfilesForEntireTraversal.getCoverageByAggregationType(partition), partition, sortedSampleLists);
    }

    /**
     * Write output summary histogram based metrics. Note that this should be invoked on a coverage partitioner that has
     * been updated for every locus across the traversal intervals.
     *
     * @param coverageProfilesForEntireTraversal DepthOfCoveragePartitionedDataStore object corresponding to the entire tool traversal.
     * @param partition                          Partition corresponding to the data to be written
     * @param sortedSamples                      sorted list of samples for this partition
     */
    public void writeCumulativeOutputSummaryFiles(final DepthOfCoveragePartitionedDataStore coverageProfilesForEntireTraversal, final DoCOutputType.Partition partition, final List<String> sortedSamples) {
        outputPerSampleCumulativeStatisticsForPartition(getCorrectOutputWriter(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.statistics), coverageProfilesForEntireTraversal.getCoverageByAggregationType(partition), sortedSamples);
        outputCumulativeSummaryForPartition(getCorrectOutputWriter(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.summary), coverageProfilesForEntireTraversal.getCoverageByAggregationType(partition), sortedSamples);
    }

    /**
     * Write out the interval summary statistics. Note that this method expects as input that the provided table has had
     * {@link CoverageUtils#updateTargetTable(int[][], DepthOfCoverageStats)} called on it exactly once for each interval
     * summarized in this traversal.
     *
     * @param partition                Partition corresponding to the data to be written
     * @param nTargetsByAvgCvgBySample Target sample coverage histogram for the given partition to be written out
     * @param binEndpoints             Bins endpoints used in the construction of of the provided histogram
     */
    public void writeOutputIntervalStatistics(final DoCOutputType.Partition partition, final int[][] nTargetsByAvgCvgBySample, final int[] binEndpoints) {
        SimpleCSVWriterWrapperWithHeader output = getCorrectOutputWriter(partition, DoCOutputType.Aggregation.interval, DoCOutputType.FileType.statistics);
        printIntervalTable(output, nTargetsByAvgCvgBySample, binEndpoints);
    }

    /**
     * Write out the gene statistics. Note that this method expects as input that the provided table has had
     * {@link CoverageUtils#updateTargetTable(int[][], DepthOfCoverageStats)} called on it exactly once for each interval
     * summarized in this traversal.
     *
     * @param nTargetsByAvgCvgBySample Target sample coverage histogram for the given partition to be written out
     * @param binEndpoints             Bins endpoints used in the construction of of the provided histogram
     */
    public void writeOutputGeneStatistics(final int[][] nTargetsByAvgCvgBySample, final int[] binEndpoints) {
        SimpleCSVWriterWrapperWithHeader output = getCorrectOutputWriter(DoCOutputType.Partition.sample, DoCOutputType.Aggregation.gene, DoCOutputType.FileType.statistics);
        printIntervalTable(output, nTargetsByAvgCvgBySample, binEndpoints);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // File writing methods
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Writes out two closely related tables, PARTITION_cumulative_coverage_counts and PARTITION_cumulative_coverage_proportions.
     * A line of PARTITION_coverage_counts stores a count of the number of loci traversed that had at least X coverage according to the
     * coverage histogram defined in the header. PARTITION_coverage_proportions corresponds to the same data except presented as a
     * percentage relative to the total number of loci traversed by the tool.
     * <p>
     * How the table is structured:
     * Number_of_sources, (for each histogram left endpoint X) +[gte_x]
     * <p>
     * Whereas a row may be recorded as follows:
     * At_least_N_samples,13,4,4,4,3,1,1,1,0,etc...
     *
     * @param countsOutput      SimpleCSVWriterWrapperWithHeader object to write the cumulative coverage counts into
     * @param proportionsOutput SimpleCSVWriterWrapperWithHeader object to write the count proportions into
     * @param stats             DepthOfCoverageStats object corresponding to the partitioning in partitionType
     * @param partitionType     partitioning this output data corresponds to
     * @param sortedSampleList  sorted list of sample names for ordering
     */
    private void outputPerLocusCumulativeSummaryAndStatistics(SimpleCSVWriterWrapperWithHeader countsOutput, SimpleCSVWriterWrapperWithHeader proportionsOutput,
                                                              DepthOfCoverageStats stats, DoCOutputType.Partition partitionType, List<String> sortedSampleList) {
        int[] endpoints = stats.getEndpoints();
        int samples = stats.getHistograms().size();

        long[][] baseCoverageCumDist = stats.getLocusCounts();

        // rows - # of samples
        // columns - depth of coverage

        List<String> headerLines = Lists.newArrayList(partitionType == DoCOutputType.Partition.readgroup ? "read_group" : partitionType.toString(), "gte_0");

        for (int d : endpoints) {
            headerLines.add("gte_" + d);
        }

        countsOutput.addHeaderLine(headerLines);
        proportionsOutput.addHeaderLine(headerLines);

        // Fill out the lines of the table to the writer based on the source tables
        for (int row = 0; row < samples; row++) {
            SimpleCSVWriterWrapperWithHeader.SimpleCSVWriterLineBuilder lineBuilder = countsOutput.getNewLineBuilder();
            lineBuilder.setColumn(0, "NSamples_" + Integer.toString( row + 1));
            for (int col = 0; col < baseCoverageCumDist[0].length; col++) {
                lineBuilder.setColumn(col + 1, Long.toString(baseCoverageCumDist[row][col]));
            }
            lineBuilder.buildAndWriteLine();
        }

        for (String sample : sortedSampleList) {
            SimpleCSVWriterWrapperWithHeader.SimpleCSVWriterLineBuilder lineBuilder = proportionsOutput.getNewLineBuilder();
            lineBuilder.setColumn(0, sample);
            double[] coverageDistribution = stats.getCoverageProportions(sample);
            for (int bin = 0; bin < coverageDistribution.length; bin++) {
                lineBuilder.setColumn(bin + 1, DOUBLE_FORMAT_2PLACES.format( coverageDistribution[bin]));
            }
            lineBuilder.buildAndWriteLine();
        }
    }

    /**
     * Constructs the header for locus summary type files. Note these output files take the form:
     * title, total_coverage, average_coverage, (for each sample) +[sample_total_cvg, sample_mean_cvg, sample_granular_Q1, sample_granular_median, sample_granular_Q3]
     *
     * @param title      Locus type to store
     * @param allSamples Samples associated with this partition in sorted order
     * @return String corresponding to the interval summary file columns in order
     */
    private List<String> getIntervalSummaryHeader(final String title, final Collection<String> allSamples) {
        List<String> summaryHeader = Lists.newArrayList(title, "total_coverage", "average_coverage");

        for (String sample : allSamples) {
            for (String suffix : summaryHeaderSampleSuffixes) {
                summaryHeader.add(sample + suffix);
            }
        }
        return summaryHeader;
    }

    // Sample suffixes that correspond to columns in the locus summary header
    private List<String> getSampleSuffixes(final List<Integer> coverageThresholds) {
        final List<String> tmp = Lists.newArrayList("_total_cvg", "_mean_cvg", "_granular_Q1", "_granular_median", "_granular_Q3");
        for (int thresh : coverageThresholds) {
            tmp.add("_%_above_" + thresh);
        }
        return Collections.unmodifiableList(tmp);
    }


    /**
     * Constructs the interval (typically gene) summary line, which consists of:
     * - A locus and the total/average coverage across all loci covered by that interval
     * - The total/average coverage for each sample in the partition at that locus
     * - The median, 1st, and 3rd quartile locus coverage statistics for each sample
     */
    private void printIntervalSummaryLine(final SimpleCSVWriterWrapperWithHeader outputWriter, final String locusName, final DepthOfCoverageStats stats, final List<String> sortedSamples) {
        SimpleCSVWriterWrapperWithHeader.SimpleCSVWriterLineBuilder lineBuilder = outputWriter.getNewLineBuilder();
        int[] bins = stats.getEndpoints();

        lineBuilder.setColumn(0, locusName)
                .setColumn("total_coverage", Long.toString(stats.getTotalCoverage()))
                .setColumn("average_coverage", DOUBLE_FORMAT_2PLACES.format(stats.getTotalMeanCoverage()));

        // each sample is in the order +[sample_total_cvg, sample_mean_cvg, sample_granular_Q1, sample_granular_median, sample_granular_Q3] so we set colums accordingly
        for (String s : sortedSamples) {
            int sIdx = outputWriter.getIndexForColumn(s + "_total_cvg");
            int median = CoverageUtils.getQuantile(stats.getHistograms().get(s), 0.5);
            int q1 = CoverageUtils.getQuantile(stats.getHistograms().get(s), 0.25);
            int q3 = CoverageUtils.getQuantile(stats.getHistograms().get(s), 0.75);
            lineBuilder.setColumn(sIdx, Long.toString(stats.getTotals().get(s)))
                    .setColumn(sIdx + 1, DOUBLE_FORMAT_2PLACES.format(stats.getMeans().get(s)))
                    .setColumn(sIdx + 2, formatBin(bins, q1))
                    .setColumn(sIdx + 3, formatBin(bins, median))
                    .setColumn(sIdx + 4, formatBin(bins, q3));

            for (int thresh : coverageThresholds) {
                lineBuilder.setColumn(s + "_%_above_" + thresh, DOUBLE_FORMAT_1PLACE.format(CoverageUtils.getPctBasesAbove(stats.getHistograms().get(s), stats.value2bin(thresh))));
            }
        }
        lineBuilder.buildAndWriteLine();
    }

    /**
     * Writes the "PARTITION_statistics" file. This file is a histogram of the total coverage in each sample at each locus.
     * <p>
     * The lines are organized as follows:
     * Source_of_reads, (for each histogram bin) +[from_X_to_Y]
     *
     * @param output        Output writer to write into
     * @param stats         DepthOfCoverageStats object constructed for this partition
     * @param sortedSamples Sorted list of all samples, to be used for ordering the output
     */
    private void outputPerSampleCumulativeStatisticsForPartition(final SimpleCSVWriterWrapperWithHeader output, final DepthOfCoverageStats stats, final List<String> sortedSamples) {
        int[] leftEnds = stats.getEndpoints();

        // Construct the header for the histogram output
        List<String> headerColumns = Lists.newArrayList("Source_of_reads");

        headerColumns.add("from_0_to_" + leftEnds[0] + ")");
        for (int i = 1; i < leftEnds.length; i++) {
            headerColumns.add("from_" + leftEnds[i - 1] + "_to_" + leftEnds[i] + ")");
        }
        headerColumns.add("from_" + leftEnds[leftEnds.length - 1] + "_to_inf");
        output.addHeaderLine(headerColumns);

        // Print out the histogram for every line
        Map<String, long[]> histograms = stats.getHistograms();
        for (String sample : sortedSamples) {
            SimpleCSVWriterWrapperWithHeader.SimpleCSVWriterLineBuilder lineBuilder = output.getNewLineBuilder();
            lineBuilder.setColumn("Source_of_reads", "sample_" + sample);
            long[] histForSample = histograms.get(sample);
            for (int i = 0; i < histForSample.length; i++) {
                lineBuilder.setColumn(i + 1, Long.toString(histForSample[i])); // +1 on the index here because the histogram starts at index 1
            }
            lineBuilder.buildAndWriteLine();
        }
    }

    /**
     * @param output        Output writer to write into
     * @param stats         DepthOfCoverageStats object constructed for this partition
     * @param sortedSamples Sorted list of all samples, to be used for ordering the output
     */
    private void outputCumulativeSummaryForPartition(final SimpleCSVWriterWrapperWithHeader output, final DepthOfCoverageStats stats, final List<String> sortedSamples) {
        List<String> headerLines = Lists.newArrayList(CUMULATIVE_SUMMARY_OUTPUT_LINES);

        for (int thresh : coverageThresholds) {
            headerLines.add("%_bases_above_" + thresh);
        }
        output.addHeaderLine(headerLines);

        Map<String, long[]> histograms = stats.getHistograms();
        Map<String, Double> means = stats.getMeans();
        Map<String, Long> totals = stats.getTotals();
        int[] leftEnds = stats.getEndpoints();

        // Write out a line for each sample
        for (String sample : sortedSamples) {
            SimpleCSVWriterWrapperWithHeader.SimpleCSVWriterLineBuilder lineBuilder = output.getNewLineBuilder();

            long[] histogram = histograms.get(sample);
            int median = CoverageUtils.getQuantile(histogram, 0.5);
            int q1 = CoverageUtils.getQuantile(histogram, 0.25);
            int q3 = CoverageUtils.getQuantile(histogram, 0.75);
            // if any of these are larger than the higest bin, put the median as in the largest bin
            median = median == histogram.length - 1 ? histogram.length - 2 : median;
            q1 = q1 == histogram.length - 1 ? histogram.length - 2 : q1;
            q3 = q3 == histogram.length - 1 ? histogram.length - 2 : q3;

            lineBuilder.setColumn(0, sample)
                    .setColumn(1, Long.toString(totals.get(sample)))
                    .setColumn(2, DOUBLE_FORMAT_2PLACES.format(means.get(sample)))
                    .setColumn(3, Integer.toString(leftEnds[q3]))
                    .setColumn(4, Integer.toString(leftEnds[median]))
                    .setColumn(5, Integer.toString(leftEnds[q1]));

            for (int thresh : coverageThresholds) {
                lineBuilder.setColumn("%_bases_above_" + thresh, DOUBLE_FORMAT_1PLACE.format(CoverageUtils.getPctBasesAbove(histogram, stats.value2bin(thresh))));
            }

            lineBuilder.buildAndWriteLine();
        }

        // Write out one final line corresponding to the total coverage
        SimpleCSVWriterWrapperWithHeader.SimpleCSVWriterLineBuilder lineBuilder = output.getNewLineBuilder();
        lineBuilder.setColumn(0, "Total")
                .setColumn(1, Long.toString(stats.getTotalCoverage()))
                .setColumn(2, DOUBLE_FORMAT_2PLACES.format(stats.getTotalMeanCoverage()))
                .fill("N/A").buildAndWriteLine();
    }

    /**
     * Produces the "PARTITION_interval_statistics" table, in which each row of the table corresponds to a count of the number of samples,
     * and each column tracks how many of the specified intervals had at least the given histogram value median depth across their spans.
     * <p>
     * How the table is structured:
     * Number_of_sources, (for each histogram left endpoint X) +[depth>=x]
     * for example:
     * Number_of_sources, depth>=0, depth>=1, depth>=2, depth>=3, depth>=4, depth>=5, depth>=6...
     * <p>
     * Whereas a row may be recorded as follows:
     * At_least_N_samples,13,4,4,4,3,1,1,1,0,etc...
     * Assuming this example corresponds to the example header, this row can be interpreted to mean that across all
     * the interval targets provided to the tool, 4 of them had a median depth of at least one base in N or more of the
     * samples in this partition.
     * <p>
     * See {@link CoverageUtils#updateTargetTable(int[][], DepthOfCoverageStats)} for more details as to how the table is constructed.
     *
     * @param output        Writer to output the table into
     * @param intervalTable Interval table object to summarize
     * @param cutoffs       histogram bin left endpoints
     */
    private void printIntervalTable(SimpleCSVWriterWrapperWithHeader output, int[][] intervalTable, int[] cutoffs) {
        List<String> columns = Lists.newArrayList("Number_of_sources", "depth>=0");
        for (int col = 0; col < intervalTable[0].length - 1; col++) {
            columns.add("depth>=" + cutoffs[col]);
        }
        output.addHeaderLine(columns);

        // Fill out the lines of the table to the writer based on the source tables
        for (int row = 0; row < intervalTable.length; row++) {
            SimpleCSVWriterWrapperWithHeader.SimpleCSVWriterLineBuilder lineBuilder = output.getNewLineBuilder();
            lineBuilder.setColumn(0, "At_least_" + Integer.toString( row + 1) + "_samples");
            for (int col = 0; col < intervalTable[0].length; col++) {
                lineBuilder.setColumn(col + 1, Integer.toString(intervalTable[row][col]));
            }
            lineBuilder.buildAndWriteLine();
        }
    }

    private String formatBin(int[] bins, int quartile) {
        if (quartile >= bins.length) {
            return ">" + Integer.toString(bins[bins.length - 1]);
        } else if (quartile < 0) {
            return "<" + Integer.toString(bins[0]);
        } else {
            return Integer.toString(bins[quartile]);
        }
    }

    /**
     * @param counts           per-nucleotide base counts for a locus, will fill in 0s if the list is null due to empty samples/sites
     * @param includeDeletions whether to include "D:##" in the output string
     * @return writer string for base counts
     */
    private String getBaseCountsString(int[] counts, boolean includeDeletions) {
        if (counts == null) {
            counts = new int[6];
        }
        StringBuilder s = new StringBuilder();
        int nbases = 0;
        for (byte b : BaseUtils.BASES_EXTENDED) {
            nbases++;
            if (includeDeletions || b != BaseUtils.Base.D.base) {
                s.append((char) b);
                s.append(":");
                s.append(counts[BaseUtils.extendedBaseToBaseIndex(b)]);
                if (nbases < 6) {
                    s.append(" ");
                }
            }
        }

        return s.toString();
    }

    @Override
    public void close() {
        try {
            for (SimpleCSVWriterWrapperWithHeader stream : outputs.values()) {
                stream.close();
            }
        } catch (IOException e) {
            throw new GATKException("Error closing output files:", e);
        }
    }


    /**
     * A simple helper class wrapper around CSVWriter that has the ingrained concept of a header line with indexed fields
     * <p>
     * Why didn't I use a tableWriter here? Who really holds the patent on the wheel anyway?
     */
    static class SimpleCSVWriterWrapperWithHeader implements Closeable {
        private int expectedColumns;
        private Map<String, Integer> headerMap = null;
        private CSVWriter outputWriter;

        // The current incomplete line in the writer.
        private SimpleCSVWriterLineBuilder currentLine = null;

        private SimpleCSVWriterWrapperWithHeader(Writer writer, char separator) {
            outputWriter = new CSVWriter(writer, separator);
        }

        /**
         * Provides a header line to the CSV output file. Note that this will throw an exception if all header lines
         * are not unique as it attempts to create an index for the provided header lines for convenience when building
         * rows of the CSV.
         *
         * @param columns Ordered list of header lines to be built into the CSV
         */
        private void addHeaderLine(final List<String> columns) {
            if (headerMap != null) {
                throw new GATKException("Should not be adding multiple header lines to a file");
            }
            outputWriter.writeNext(columns.toArray(new String[0]), false);
            expectedColumns = columns.size();

            // Create the mapping between header and column
            headerMap = new HashMap<>();
            for (int i = 0; i < columns.size(); i++) {
                final String columnHeading = columns.get(i);
                Utils.nonNull(columnHeading);
                if (headerMap.putIfAbsent(columnHeading, i) != null)  {
                    throw new GATKException("Only allow unique column headings");
                }
            }
        }

        private void writeLine(String[] line) {
            outputWriter.writeNext(line, false);
            currentLine = null;
        }

        /**
         * Builds a new SimpleCSVWriterLineBuilder and writes out the previous line if it exists.
         *
         * @return a blank SimpleCSVWriterLineBuilder to allow for defining the next line
         */
        private SimpleCSVWriterLineBuilder getNewLineBuilder() {
            if (headerMap == null) {
                throw new GATKException("Cannot construct line without first setting the header line");
            }
            if (currentLine != null) {
                currentLine.buildAndWriteLine();
            }
            currentLine = new SimpleCSVWriterLineBuilder(this, expectedColumns);
            return currentLine;
        }

        /**
         * @param column header line to get index for
         * @return zero based index corresponding to that header string, throws an exception if the headerline doesn't exist
         */
        private int getIndexForColumn(String column) {
            Utils.nonNull(headerMap, "Cannot request column index if the header has not been specified");
            int index = headerMap.get(column);
            Utils.nonNull(index, "Requested column " + column + " does not exist in the provided header");
            return index;
        }

        @Override
        public void close() throws IOException {
            if (currentLine != null) {
                currentLine.buildAndWriteLine();
            }
            outputWriter.close();
        }

        /**
         * Helper to allow for incremental construction of a body line using either indexes or column headings
         * <p>
         * Calling buildAndWriteLine() will cause the line to be written out into the underlying CSV writer in its current state. Doing
         * so will result in a validation call where an exception will be thrown if any columns of the current line have
         * not been defined. fill() can be used to provide a default value for undefined columns.
         */
        private static class SimpleCSVWriterLineBuilder {
            SimpleCSVWriterWrapperWithHeader thisBuilder;
            String[] lineToBuild;
            boolean hasBuilt = false;

            SimpleCSVWriterLineBuilder(SimpleCSVWriterWrapperWithHeader me, int lineLength) {
                thisBuilder = me;
                lineToBuild = new String[lineLength];
            }

            /**
             * @param index Column index to be set
             * @param value Value to be placed into the line
             */
            private SimpleCSVWriterLineBuilder setColumn(final int index, final String value) {
                checkAlteration();
                lineToBuild[index] = value;
                return this;
            }

            /**
             * @param heading Column heading to be set
             * @param value   Value to be placed into the line
             */
            private SimpleCSVWriterLineBuilder setColumn(final String heading, final String value) {
                int index = thisBuilder.getIndexForColumn(heading);
                return setColumn(index, value);
            }

            /**
             * Fills in every empty column of the pending line with the provided value
             */
            private SimpleCSVWriterLineBuilder fill(final String filling) {
                checkAlteration();
                for (int i = 0; i < lineToBuild.length; i++) {
                    if (lineToBuild[i] == null) {
                        lineToBuild[i] = filling;
                    }
                }
                return this;
            }

            /**
             * Constructs the line and writes it out to the output
             */
            private void buildAndWriteLine() {
                Utils.validate(!Arrays.stream(lineToBuild).anyMatch(Objects::isNull), "Attempted to construct an incomplete line, make sure all columns are filled");
                thisBuilder.writeLine(lineToBuild);
                hasBuilt = true;
            }

            // Throw an exception if we try to alter an already written out line
            private void checkAlteration() {
                Utils.validate(!hasBuilt, "Cannot make alterations to an already written out CSV line");
            }
        }
    }
}