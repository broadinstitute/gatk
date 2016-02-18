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

    @Argument(doc = "Should an output plot be created")
    public boolean PRODUCE_FILES = false;

    @Argument(doc = "A local path to file to write metrics to;",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    public String OUTPUT;

    @Argument(doc="PDF File to write insert size Histogram chart to.",
            shortName="H")
    public File HISTOGRAM_FILE;

    @Argument(doc="Similar to Picard CollectInsertSizeMetrics::DEVIATIONS")
    public double DEVIATIONS_TOL = 10;

    @Argument(doc="The level(s) at which to accumulate metrics",
            shortName="LEVEL")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = EnumSet.of(MetricAccumulationLevel.ALL_READS);

    // read filtering criteria
    @Argument(shortName = "F",
            fullName = "pfReadsOnly",
            doc = "If set to true, include PF reads only.")
    public boolean pfReadsOnly = false;

    @Argument(shortName = "U",
            fullName = "mappedReadsOnly",
            doc = "If set to true, include mapped reads only.")
    public boolean mappedReadsOnly = false;

    @Argument(shortName = "Um",
            fullName = "mateMappedReadsOnly",
            doc = "If set to true, include mate-mapped reads only.")
    public boolean mateMappedReadsOnly = false;

    @Argument(shortName = "D",
            fullName = "nonDupReadsOnly",
            doc = "If set to true, inlcude non-duplicated reads only.")
    public boolean nonDupReadsOnly = false;

    @Argument(shortName = "S",
            fullName = "nonSecondaryAlignmentsOnly",
            doc = "If set to true, filter out secondary aglignments reads.")
    public boolean nonSecondaryAlignmentsOnly = false;

    @Argument(shortName = "SS",
            fullName = "nonSupplementaryAlignmentsOnly",
            doc = "If set to true, filter out supplementary aglignments reads.")
    public boolean nonSupplementaryAlignmentsOnly = false;

    @Argument(shortName = "PP",
            fullName = "properlyPariedOnly",
            doc = "If set to true, filter out supplementary aglignments reads.")
    public boolean properlyPairedOnly = false;

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
    protected void runTool(final JavaSparkContext ctx)
    {
        // two step filter: first is general, makeReadFilter, implicitly called in getReads;
        //   second is MetricsReadFilter, but provides limited capability now.
        //   So makeReadFilter is hacked.
        //   TODO: improve MetricsReadFilter that uses SAM/BAM/CRAM flags.
        final JavaRDD<GATKRead> reads = getReads();
        // final MetricsReadFilter metricsFilter = new MetricsReadFilter(pfReadsOnly, mappedReadsOnly);
        // final JavaRDD<GATKRead> filteredReads = reads.filter(read -> metricsFilter.test(read));

        collector = new ISZMetricsCollectorSpark(reads, METRIC_ACCUMULATION_LEVEL, DEVIATIONS_TOL);

        if(PRODUCE_FILES)
            saveResults(getHeaderForReads(), getReadSourceName());
    }

    @Override
    public ReadFilter makeReadFilter()
    {
        //return ReadFilterLibrary.ALLOW_ALL_READS;
        return new CustomReadFilter();
    }

    /**
     * Implements serializable read filters based on cmd line arguments provided.
     */
    private final class CustomReadFilter implements ReadFilter {
        private static final long serialVersionUID = 1L;

        @Override
        public boolean test(GATKRead read)
        {
            boolean result = false;
            if(pfReadsOnly)
                result &= !read.failsVendorQualityCheck();
            if(mappedReadsOnly)
                result &= !read.isUnmapped();
            if(mateMappedReadsOnly)
                result &= !read.mateIsUnmapped();
            if(nonDupReadsOnly)
                result &= !read.isDuplicate();
            if(nonSecondaryAlignmentsOnly)
                result &= !read.isSecondaryAlignment();
            if(nonSupplementaryAlignmentsOnly)
                result &= !read.isSupplementaryAlignment();
            if(properlyPairedOnly)
                result &= read.isProperlyPaired();

            if(MQPassingThreshold > 0.0)
                result &= (read.getMappingQuality() >= MQPassingThreshold);

            // TODO: efficiently checking mate mapping quality
            // TODO: pick only isFirstOfPair/isSecondOfPair.
            // TODO: pick only isReverseStrand/mateIsReverseStrand
            // TODO: filter or not based on length==0
            return result;
        }
    }

    /**
     * Produce metrics files and histogram chart upon request
     */
    protected void saveResults(final SAMFileHeader readsHeader, final String inputFileName)
    {
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