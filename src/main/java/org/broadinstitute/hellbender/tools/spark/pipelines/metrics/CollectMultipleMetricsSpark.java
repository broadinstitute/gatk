package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.metrics.Header;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.metrics.MetricsArgumentCollection;
import org.broadinstitute.hellbender.metrics.QualityYieldMetricsArgumentCollection;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.util.*;

/**
 * Tool that instantiates and executes multiple metrics programs using a single RDD.
 */
@CommandLineProgramProperties(
        summary = "Takes an input SAM/BAM/CRAM file and reference sequence and runs one or more " +
                "metrics modules at the same time to cut down on I/O. Currently all programs are run with " +
                "default options and fixed output extensions, but this may become more flexible in future.",
        oneLineSummary = "A \"meta-metrics\" calculating program that produces multiple metrics for the provided SAM/BAM/CRAM file",
        programGroup = SparkProgramGroup.class
)
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

    @Argument(shortName="LEVEL", doc="The level(s) at which to accumulate metrics. ", optional = true)
    public Set<MetricAccumulationLevel> metricAccumulationLevel = EnumSet.of(MetricAccumulationLevel.ALL_READS);

    @Argument(doc = "List of metrics collectors to apply during the pass through the SAM file.")
    public List<Collectors> collectors = new ArrayList<>();

    private interface CollectorProvider {
        /**
         * For each collector type, provide a type-safe, collector-specific method
         * that creates and populates an instance of the  collector's input argument
         * class; creates an instance of the collector; initializes the collector
         * with the arguments, and returns the initialized collector.
         */
        MetricsCollectorSpark<? extends MetricsArgumentCollection> createCollector(
                final String outputBaseName,
                final Set<MetricAccumulationLevel> metricAccumulationLevel,
                final List<Header> defaultHeaders
        );
    }

    // Enum of Collector types that CollectMultipleMetrics can support.
    public static enum Collectors implements CollectorProvider {
        CollectInsertSizeMetrics {
            @Override
            public MetricsCollectorSpark<? extends MetricsArgumentCollection> createCollector(
                    final String outputBaseName,
                    final Set<MetricAccumulationLevel> metricAccumulationLevel,
                    final List<Header> defaultHeaders)
            {
                // disambiguate this collector's output files from the other collectors
                final String localBaseName = outputBaseName + ".insertMetrics";

                final InsertSizeMetricsArgumentCollection isArgs = new InsertSizeMetricsArgumentCollection();
                isArgs.output = new File(localBaseName + ".txt");
                isArgs.histogramPlotFile = localBaseName + ".pdf";
                isArgs.useEnd = InsertSizeMetricsArgumentCollection.EndToUse.SECOND;
                isArgs.metricAccumulationLevel = metricAccumulationLevel;

                final InsertSizeMetricsCollectorSpark collector = new InsertSizeMetricsCollectorSpark();
                collector.initialize(isArgs, defaultHeaders);

                return collector;
            }
        },
        CollectQualityYieldMetrics {
            @Override
            public MetricsCollectorSpark<? extends MetricsArgumentCollection> createCollector(
                    final String outputBaseName,
                    final Set<MetricAccumulationLevel> metricAccumulationLevel,
                    final List<Header> defaultHeaders)
            {
                // disambiguate this collector's output files from the other collectors
                final String localBaseName = outputBaseName + ".qualityYieldMetrics";

                final QualityYieldMetricsArgumentCollection qyArgs = new QualityYieldMetricsArgumentCollection();
                qyArgs.output = new File(localBaseName + ".txt");

                final QualityYieldMetricsCollectorSpark collector = new QualityYieldMetricsCollectorSpark();
                collector.initialize(qyArgs, defaultHeaders);

                return collector;
            }
        }
    }

    @Override
    protected void runTool( final JavaSparkContext ctx ) {
        final JavaRDD<GATKRead> unFilteredReads = getUnfilteredReads();

        if (collectors.isEmpty()) { // run all collectors
            collectors.addAll(Arrays.asList(Collectors.values()));
        }
        if (collectors.size() > 1) {
            // if there is more than one collector to run, cache the
            // unfiltered RDD so we don't recompute it
            unFilteredReads.cache();
        }

        for (final CollectorProvider provider : collectors) {
            MetricsCollectorSpark<? extends MetricsArgumentCollection> metricsCollector =
                    provider.createCollector(
                        outputBaseName,
                        metricAccumulationLevel,
                        getDefaultHeaders()
                    );

            validateCollector(metricsCollector, collectors.get(collectors.indexOf(provider)).name());

            // Execute the collector's lifecycle
            ReadFilter readFilter =
                    metricsCollector.getReadFilter(getHeaderForReads());
            metricsCollector.collectMetrics(
                    unFilteredReads.filter(r -> readFilter.test(r)),
                    getHeaderForReads(),
                    getReadSourceName(),
                    getAuthHolder()
            );
            metricsCollector.finishCollection();
        }
    }

    // Validate requiresReference and expected sort order
    private void validateCollector(final MetricsCollectorSpark<?> collector, final String name) {
        if (collector.requiresReference()) {
            throw new UnsupportedOperationException("Requires reference for collector not yet implemented");
        }
        ReadUtils.validateExpectedSortOrder(
                getHeaderForReads().getSortOrder(),
                collector.getExpectedSortOrder(),
                false,
                name);
    }

}
