package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.filters.AlignmentAgreesWithHeaderReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.metrics.MetricsUtils;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetrics;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.util.*;

@CommandLineProgramProperties(
        summary        = "Program to collect insert size distribution information in SAM/BAM file(s)",
        oneLineSummary = "CollectInsertSizeDistr on Spark",
        programGroup = SparkProgramGroup.class)
public final class CollectISizeDistributionSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "Should an output plot be created?")
    public boolean PRODUCE_FILES = false;

    @Argument(doc = "A local path to file to write insert size metrics to.",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    public String OUTPUT;

    @Argument(doc="PDF File to write histogram chart to.",
            shortName="H")
    public File HISTOGRAM_FILE;

    @Argument(doc="Similar to Picard CollectInsertSizeMetrics::DEVIATIONS.")
    public double DEVIATIONS_TOL = 10;

    @Argument(doc="The level(s) at which to accumulate metrics.",
            shortName="LEVEL")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = EnumSet.of(MetricAccumulationLevel.ALL_READS);

    // read filtering criteria

    @Argument(shortName = "PP",
            fullName = "properlyPariedOnly",
            doc = "If set to true, filter out pairs of reads that are not properly oriented. Default value: true.")
    public boolean properlyPairedOnly = true;

    @Argument(shortName = "U",
            fullName = "mappedReadsOnly",
            doc = "If set to true, include mapped reads only. Default value: true.")
    public boolean mappedReadsOnly = true;

    @Argument(shortName = "Um",
            fullName = "mateMappedReadsOnly",
            doc = "If set to true, include mate-mapped reads only. Default value: true.")
    public boolean mateMappedReadsOnly = true;

    @Argument(shortName = "D",
            fullName = "nonDupReadsOnly",
            doc = "If set to true, include non-duplicated reads only. Default value: true.")
    public boolean nonDupReadsOnly = true;

    @Argument(shortName = "S",
            fullName = "nonSecondaryAlignmentsOnly",
            doc = "If set to true, filter out secondary alignments reads. Default value: true.")
    public boolean nonSecondaryAlignmentsOnly = true;

    @Argument(shortName = "SS",
            fullName = "nonSupplementaryAlignmentsOnly",
            doc = "If set to true, filter out supplementary alignments reads. Default value: true.")
    public boolean nonSupplementaryAlignmentsOnly = true;

    @Argument(shortName = "F",
            fullName = "whichEndOfPairToUse",
            doc = "Which end of pairs to use for collecting information. Possible values:{0,1,2}, where 1 and 2 " +
                    "stands for first or second end respectively. Option 0 picks up information from 1st end " +
                    "when both ends are available, and pick either one when only one end is available. " +
                    "(Remember, in SV analysis, a truncated region may be investigated, so this is possible.)" +
                    " Unused if only paired reads are used (properlyPairedOnly==true).")
    public int whichEndOfPairToUse = 1;

    @Argument(shortName = "Q",
            fullName = "MQPassingThreshold",
            doc = "If set non-zero value, only include reads passing certain mapping quality threshold.")
    public double MQPassingThreshold = 0.0;

    @Argument(shortName = "Qm",
            fullName = "mateMQPassingThreshold",
            doc = "If set to non-zero value, only include reads with mate passing certain mapping quality threshold.")
    public double mateMQPassingThreshold = 0.0; // TODO: efficient way to test


    // path to Picard R script for producing histograms in PDF files.
    private static final String R_SCRIPT = "insertSizeHistogram.R";

    // should this be static?
    // Calculates InsertSizeMetrics for all METRIC_ACCUMULATION_LEVELs provided
    private static ISZMetricsCollectorSpark collector;

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaRDD<GATKRead> reads = getReads();

        collector = new ISZMetricsCollectorSpark(reads, METRIC_ACCUMULATION_LEVEL, getHeaderForReads().getReadGroups(), DEVIATIONS_TOL);

        if(PRODUCE_FILES)
            saveResults(getHeaderForReads(), getReadSourceName());
    }

    @Override
    public ReadFilter makeReadFilter() { // return type serializable. implicitly called in getReads(): behavior from GATKTools
        final boolean filters[] = {properlyPairedOnly, mappedReadsOnly, mateMappedReadsOnly,
                                   nonDupReadsOnly, nonSecondaryAlignmentsOnly, nonSupplementaryAlignmentsOnly};

        final SVCustomReadFilter svfilter = new SVCustomReadFilter(filters, whichEndOfPairToUse,
                                                                   MQPassingThreshold, mateMQPassingThreshold);

        final AlignmentAgreesWithHeaderReadFilter alignmentAgreesWithHeader = new AlignmentAgreesWithHeaderReadFilter(getHeaderForReads());

        return svfilter.and(alignmentAgreesWithHeader); // Is this bad practice (type mixture)?
    }

    /**
     * Produce metrics files and histogram chart upon request
     */
    protected void saveResults(final SAMFileHeader readsHeader, final String inputFileName) {
        // write text-based metrics
        if(OUTPUT != null){
            final MetricsFile<InsertSizeMetrics, Integer> metrics = getMetricsFile();

            if (metrics.getAllHistograms().isEmpty()) { // early return if nothing to produce
                logger.warn("No valid reads found in input file.");
                return;
            }

            collector.produceMetricsFile(metrics);

            MetricsUtils.saveMetrics(metrics, OUTPUT, getAuthHolder());
        }

        // Now run Picard R script to generate graphics
        if(HISTOGRAM_FILE != null){
            IOUtil.assertFileIsWritable(HISTOGRAM_FILE);

            /* Picard R script for this purpose doesn't produce subtitle based on library
            // If we're working with a single library, assign that library's name as a suffix to the plot title
            final List<SAMReadGroupRecord> readGroups = readsHeader.getReadGroups();

            //A subtitle for the plot, usually corresponding to a library.
            String plotSubtitle = "";
            if (readGroups.size() == 1) {
                plotSubtitle = StringUtil.asEmptyIfNull(readGroups.get(0).getLibrary());
            }*/

            final RScriptExecutor executor = new RScriptExecutor();
            executor.addScript(new Resource(CollectISizeDistributionSpark.R_SCRIPT, ISZMetricsCollectorSpark.class));
            executor.addArgs(OUTPUT,                            // text-based metrics file
                             HISTOGRAM_FILE.getAbsolutePath(),  // PDF graphics file
                             inputFileName);                    // input bam file
                             //, plotSubtitle); // Picard R script for this program doesn't take subtitle argument
            executor.exec();
        }
    }
}