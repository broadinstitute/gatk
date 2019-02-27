package org.broadinstitute.hellbender.tools.examples.metrics.multi;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.metrics.Header;

import org.apache.spark.api.java.JavaRDD;

import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.spark.pipelines.metrics.MetricsCollectorSparkTool;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;

/**
 * Example Spark tool for collecting multi-level metrics.
 *
 * Collector tools consist of a thin facade that delegates to an embedded instance of a shared collector
 * implementation (in this case {@link ExampleMultiMetricsCollectorSpark})that can also be used in other
 * contexts. This example illustrates how to use the shared implementation in a standalone Spark tool.
 * The embedded collector implementation could also be used in other contexts,
 * such as from within {@link org.broadinstitute.hellbender.tools.spark.pipelines.metrics.CollectMultipleMetricsSpark}
 *
 * See {@link org.broadinstitute.hellbender.tools.spark.pipelines.metrics.MetricsCollectorSpark},
 * {@link org.broadinstitute.hellbender.tools.spark.pipelines.metrics.MetricsCollectorSparkTool},
 * {@link org.broadinstitute.hellbender.metrics.MultiLevelCollector}, and
 * <a href="http://github.com/broadinstitute/gatk/wiki/GATK-Metrics-Collectors">GATK Collectors/a>
 * for more information about GATK collectors.
 */
@CommandLineProgramProperties(
        summary        = "Program to collect example multi-level metrics in SAM/BAM/CRAM file(s)",
        oneLineSummary = "Collect example multi-level metrics on Spark",
        programGroup   = ExampleProgramGroup.class,
        omitFromCommandLine = true)
public final class ExampleCollectMultiMetricsSpark
        extends MetricsCollectorSparkTool<ExampleMultiMetricsArgumentCollection> {

    private static final long serialVersionUID = 1L;

    /**
     * Embedded instance of the collector's argument collection.
     */
    @ArgumentCollection
    private ExampleMultiMetricsArgumentCollection exampleArgs = new ExampleMultiMetricsArgumentCollection();

    /**
     * Embedded instance of the actual collector to which we will delegate
     */
    private ExampleMultiMetricsCollectorSpark exampleMultiCollector = new ExampleMultiMetricsCollectorSpark();

    /**
     * Return the embedded argument collection. Called by the base class; the results are passed back
     * to the {@link #initialize(ExampleMultiMetricsArgumentCollection, SAMFileHeader, List)} method.
     *
     * Note: In a standalone tool context such as this, the argument values are populated by the user
     * via the command line, and are maintained by the tool. In other contexts, the arguments that
     * are passed to the actual collector could be obtained from other sources. This method exists to
     * allow the tool framework to use the same collector initialization protocol for all contexts,
     * though in the tool context it represents the degenerate case since the arguments are first
     * obtained from the tool, and then immediately passed back.
     *
     * @return ExampleMultiMetricsArgumentCollection
     */
    @Override
    protected ExampleMultiMetricsArgumentCollection getInputArguments() {
        return exampleArgs;
    }

    /**
     * Return the sort order required/expected by this collector.
     * @return SAMFileHeader.SortOrder for this collector
     */
    @Override
    protected SortOrder getExpectedSortOrder() { return exampleMultiCollector.getExpectedSortOrder(); }

    /**
     * Initialize the collector with it's input arguments.
     *
     * @param inputArgs The input arguments for this collector. The {@link MetricsCollectorSparkTool} base class
     *                  will have previously obtained these through a call to {@link #getInputArguments},
     *                  and passes them back here to be forwarded to the embedded collector.
     * @param samHeader SAMFileHeader for the input
     * @param defaultHeaders Default headers for this tool to optionally be used in any resulting metrics file.
     */
    @Override
    protected void initialize(
            final ExampleMultiMetricsArgumentCollection inputArgs,
            final SAMFileHeader samHeader,
            final List<Header> defaultHeaders)
    {
        exampleMultiCollector.initialize(inputArgs, samHeader, defaultHeaders);
    }

    /**
     * Return the read filters required for this collector
     */
    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return exampleMultiCollector.getDefaultReadFilters();
    }

    /**
     * Execute the actual metrics collection.
     * @param filteredReads Input reads, already filtered using the filter returned by {@link ##getDefaultReadFilters(SAMFileHeader)}
     * @param samHeader SAMFileHeader for the input
     */
    @Override
    protected void collectMetrics(
            final JavaRDD<GATKRead> filteredReads,
            final SAMFileHeader samHeader) {
        exampleMultiCollector.collectMetrics(filteredReads, samHeader);
    }

    /**
     * Finish the collection process and save any results.
     * @param inputName The name of the input source for optional inclusion in the saved metrics file.
     */
    @Override
    protected void saveMetrics(final String inputName) {
        exampleMultiCollector.saveMetrics(inputName);
    }
}
