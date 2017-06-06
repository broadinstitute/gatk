package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.Header;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MetricAccumulationLevelArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.metrics.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;

/**
 * Tool that instantiates and executes multiple metrics programs using a single RDD.
 */
@CommandLineProgramProperties(
        summary = "Takes an input SAM/BAM/CRAM file and reference sequence and runs one or more Spark " +
                "metrics modules at the same time to cut down on I/O. Currently all programs are run with " +
                "default options and fixed output extensions, but this may become more flexible in future.",
        oneLineSummary = "A \"meta-metrics\" calculating program that produces multiple metrics for the provided SAM/BAM/CRAM file",
        programGroup = SparkProgramGroup.class
)
@DocumentedFeature
public final class CollectMultipleMetricsSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Argument(fullName=StandardArgumentDefinitions.ASSUME_SORTED_LONG_NAME,
                shortName = StandardArgumentDefinitions.ASSUME_SORTED_SHORT_NAME,
                doc = "If true (default), then the sort order in the header file will be ignored.")
    public boolean ASSUME_SORTED = true;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
                shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
                doc = "Base name of output files.")
    public String outputBaseName;

    @ArgumentCollection
    MetricAccumulationLevelArgumentCollection metricAccumulationLevel = new MetricAccumulationLevelArgumentCollection();

    @Argument(fullName="collectors",
            doc = "List of metrics collectors to apply during the pass through the SAM file. " +
                  "If no collectors are specified than all collectors will be run", optional = true)
    public List<SparkCollectors> userCollectors = new ArrayList<>();

    public interface SparkCollectorProvider {
        /**
         * For each collector type, provide a type-safe, collector-specific method
         * that creates and populates an instance of the  collector's input argument
         * class; creates an instance of the collector; initializes the collector
         * with the arguments, and returns the initialized collector.
         */
        MetricsCollectorSpark<? extends MetricsArgumentCollection> createCollector(
                final String outputBaseName,
                final Set<MetricAccumulationLevel> metricAccumulationLevel,
                final List<Header> defaultHeaders,
                final SAMFileHeader samHeader
        );
    }

    // Enum of Collector types that CollectMultipleMetrics can support.
    public static enum SparkCollectors implements SparkCollectorProvider {
        CollectInsertSizeMetrics {
            @Override
            public MetricsCollectorSpark<? extends MetricsArgumentCollection> createCollector(
                    final String outputBaseName,
                    final Set<MetricAccumulationLevel> metricAccumulationLevel,
                    final List<Header> defaultHeaders,
                    final SAMFileHeader samHeader)
            {
                // disambiguate this collector's output files from the other collectors
                final String localBaseName = outputBaseName + "." + InsertSizeMetrics.getUniqueNameSuffix();

                final InsertSizeMetricsArgumentCollection isArgs = new InsertSizeMetricsArgumentCollection();
                isArgs.output = localBaseName + ".txt";
                isArgs.histogramPlotFile = localBaseName + ".pdf";
                isArgs.metricAccumulationLevel.accumulationLevels = metricAccumulationLevel;

                final InsertSizeMetricsCollectorSpark collector = new InsertSizeMetricsCollectorSpark();
                collector.initialize(isArgs, samHeader, defaultHeaders);

                return collector;
            }
        },
        CollectQualityYieldMetrics {
            @Override
            public MetricsCollectorSpark<? extends MetricsArgumentCollection> createCollector(
                    final String outputBaseName,
                    final Set<MetricAccumulationLevel> metricAccumulationLevel,
                    final List<Header> defaultHeaders,
                    final SAMFileHeader samHeader)
            {
                // disambiguate this collector's output files from the other collectors
                final String localBaseName = outputBaseName + "." + QualityYieldMetrics.getUniqueNameSuffix();

                final QualityYieldMetricsArgumentCollection qyArgs = new QualityYieldMetricsArgumentCollection();
                qyArgs.output = localBaseName + ".txt";

                final QualityYieldMetricsCollectorSpark collector = new QualityYieldMetricsCollectorSpark();
                collector.initialize(qyArgs, samHeader, defaultHeaders);

                return collector;
            }
        }
    }

    /**
     * List of external collectors so that an outside developer can invoke this class programmatically
     * and provide alternative collectors to run by calling setCollectorsToRun().
     */
    private List<SparkCollectorProvider> externalCollectors = null;

    /**
     * Use this method when invoking CollectMultipleMetricsSpark programmatically to run programs other
     * than the ones available via enum. This must be called before runTool().
     */
    public void setCollectorsToRun(List<SparkCollectorProvider> externalCollectorsToRun) {
        this.externalCollectors = externalCollectorsToRun;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx ) {
        final JavaRDD<GATKRead> unFilteredReads = getUnfilteredReads();
        List<SparkCollectorProvider> collectorsToRun = getCollectorsToRun();
        if (collectorsToRun.size() > 1) {
            // if there is more than one collector to run, cache the
            // unfiltered RDD so we don't recompute it
            unFilteredReads.cache();
        }
        for (final SparkCollectorProvider provider : collectorsToRun) {
            MetricsCollectorSpark<? extends MetricsArgumentCollection> metricsCollector =
                    provider.createCollector(
                        outputBaseName,
                        metricAccumulationLevel.accumulationLevels,
                        getDefaultHeaders(),
                        getHeaderForReads()
                    );
            validateCollector(metricsCollector, collectorsToRun.get(collectorsToRun.indexOf(provider)).getClass().getName());

            // Execute the collector's lifecycle

            //Bypass the framework merging of command line filters and just apply the default
            //ones specified by the collector
            ReadFilter readFilter = ReadFilter.fromList(metricsCollector.getDefaultReadFilters(), getHeaderForReads());

            metricsCollector.collectMetrics(
                    unFilteredReads.filter(r -> readFilter.test(r)),
                    getHeaderForReads()
            );
            metricsCollector.saveMetrics(getReadSourceName());
        }
    }

    /**
     * Determine which collectors to run based on commandline args and programmatically
     * specified custom collectors.
     */
    private List<SparkCollectorProvider> getCollectorsToRun() {
        if (externalCollectors == null) {
            if (userCollectors.size() == 0) { // run them all
                return Arrays.asList(SparkCollectors.values());
            }
            else { // run the user specified collectors
                List<SparkCollectorProvider> collectorsToRun = new ArrayList<>();
                collectorsToRun.addAll(userCollectors);
                return collectorsToRun;
            }
        }
        else { // run the external collectors
            return externalCollectors;
        }
    }

    // Validate requiresReference and expected sort order
    private void validateCollector(final MetricsCollectorSpark<?> collector, final String sourceName) {
        if (collector.requiresReference()) {
            throw new UnsupportedOperationException("Requires reference for collector not yet implemented");
        }
        ReadUtils.validateExpectedSortOrder(
                getHeaderForReads().getSortOrder(),
                collector.getExpectedSortOrder(),
                false,
                sourceName);
    }

}
