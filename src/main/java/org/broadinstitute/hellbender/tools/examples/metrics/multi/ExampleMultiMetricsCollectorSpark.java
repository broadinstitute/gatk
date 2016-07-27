package org.broadinstitute.hellbender.tools.examples.metrics.multi;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.pipelines.metrics.MetricsCollectorSpark;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;

/**
 * Example implementation of a multi-level Spark metrics collector. Can be embedded within a standalone tool
 * (see @link ExampleCollectMultipleMetricsSpark}, a pipeline, or
 * {@link org.broadinstitute.hellbender.tools.spark.pipelines.metrics.CollectMultipleMetricsSpark}.
 */
public class ExampleMultiMetricsCollectorSpark
        implements MetricsCollectorSpark<ExampleMultiMetricsArgumentCollection>,
        Serializable
{
    private static final long serialVersionUID = 1L;

    // This example collector has no actual parameters that are controlled by inputArguments, but
    // this field is retained to illustrate how they pass through a typical lifecycle.
    private ExampleMultiMetricsArgumentCollection inputArgs = null;

    // Container for per-partition metrics
    private ExampleMultiMetricsCollector collector = null;

    // Container for final aggregate reduced/combined metrics
    private ExampleMultiMetricsCollector reducedResultMetrics = null;

    private MetricsFile<ExampleMultiMetrics, Integer> metricsFile = null;

    /**
     * Initialize the collector with input arguments.
     */
    @Override
    public void initialize(
            final ExampleMultiMetricsArgumentCollection inputArgs,
            final SAMFileHeader samHeader,
            final List<Header> defaultHeaders)
    {
        this.inputArgs = inputArgs;
        collector = new ExampleMultiMetricsCollector(inputArgs, samHeader);
        metricsFile = new MetricsFile<>();
        if (defaultHeaders != null) {
            defaultHeaders.stream().forEach(h -> metricsFile.addHeader(h));
        }
    }

    /**
     * Return a read filter to be used for this collector.
     *
     * @param samFileHeader
     * @return ReadFilter
     */
    @Override
    public ReadFilter getReadFilter(final SAMFileHeader samFileHeader) {
        return collector.getReadFilter(samFileHeader);
    }

    /**
     * Do he actual Metrics collection.
     * @param filteredReads The reads to be analyzed for this collector. The reads will have already
     *                      been filtered by this collector's read filter.
     * @param samHeader The SAMFileHeader associated with the reads in the input RDD.
     */
    @Override
    public void collectMetrics(
            final JavaRDD<GATKRead> filteredReads,
            final SAMFileHeader samHeader)
    {
        if (filteredReads.isEmpty()) {
            throw new GATKException("No valid reads found in input file.");
        }

        // Each partition iterates through it's subset of the RDD; the results are then
        // collected and reduced into the final aggregate metrics.
        reducedResultMetrics = filteredReads
            .mapPartitions(
                it -> {
                    it.forEachRemaining(r -> collector.acceptRecord(r.convertToSAMRecord(samHeader), null));
                    return Arrays.asList(collector);
                })
            .reduce((c1, c2) -> ExampleMultiMetricsCollector.combine(c1, c2));
    }

    /**
     * Save the metrics out to a file.
     * @param inputName name of the input from which metrics were collected; unused here
     *                  but optionally used by some collectors to label output plots, etc.
     * @param authHolder authentication info. May be null.
     */
    @Override
    public void saveMetrics(final String inputName, final AuthHolder authHolder) {
        reducedResultMetrics.saveMetrics(metricsFile, authHolder);
    }

}