package org.broadinstitute.hellbender.tools.examples.metrics.single;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.metrics.MetricsUtils;
import org.broadinstitute.hellbender.tools.spark.pipelines.metrics.MetricsCollectorSpark;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.List;

/**
 * ExampleSingleMetricsCollector for Spark.
 */
public class ExampleSingleMetricsCollectorSpark
        implements MetricsCollectorSpark<ExampleSingleMetricsArgumentCollection>, Serializable {

    private static final long serialVersionUID = 1L;

    // This example collector has no actual parameters that are controlled by input arguments, but
    // this field is retained to illustrate how such arguments are passed through a typical collector
    // lifecycle.
    private ExampleSingleMetricsArgumentCollection args = null;

    private MetricsFile<ExampleSingleMetrics, Integer> metricsFile = null;

    /**
     * Initialize the collector with input arguments;
     */
    @Override
    public void initialize(
            final ExampleSingleMetricsArgumentCollection inputArgs,
            final SAMFileHeader samHeader,
            final List<Header> defaultHeaders) {
        metricsFile = new MetricsFile<ExampleSingleMetrics, Integer>();
        defaultHeaders.stream().forEach(h -> metricsFile.addHeader(h));
        this.args = inputArgs;
    }

    /**
     * Do the actual metrics collection on the provided RDD.
     * @param filteredReads The reads to be analyzed for this collector.
     * @param samHeader The SAMFileHeader associated with the reads in the input RDD.
     */
    @Override
    public void collectMetrics(
            final JavaRDD<GATKRead> filteredReads,
            final SAMFileHeader samHeader)
    {
        final ExampleSingleMetrics metrics =
                filteredReads.aggregate(new ExampleSingleMetrics(),
                        (hgp, read) -> hgp.addRead(read),
                        (hgp1, hgp2) -> hgp1.combine(hgp2))
                        .finish();

        metricsFile.addMetric(metrics);
    }

    /**
     * Write out the results of the metrics collection
     * @param inputBaseName base name of the input file
     */
    @Override
    public void saveMetrics(final String inputBaseName) {
        MetricsUtils.saveMetrics(metricsFile, args.output);
    }

}
