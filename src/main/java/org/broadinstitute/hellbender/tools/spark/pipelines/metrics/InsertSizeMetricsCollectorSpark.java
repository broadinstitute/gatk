package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.metrics.InsertSizeMetrics;
import org.broadinstitute.hellbender.metrics.InsertSizeMetricsArgumentCollection;
import org.broadinstitute.hellbender.metrics.InsertSizeMetricsCollector;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;

/**
 * Worker class to collect insert size metrics, add metrics to file, and provides
 * accessors to stats of groups of different level.
 */
public class InsertSizeMetricsCollectorSpark implements
        MetricsCollectorSpark<InsertSizeMetricsArgumentCollection>, Serializable {
    private static final long serialVersionUID = 1L;

    private InsertSizeMetricsArgumentCollection inputArgs = null;
    private InsertSizeMetricsCollector collector = new InsertSizeMetricsCollector();

    private InsertSizeMetricsCollector resultMetrics = null;
    private MetricsFile<InsertSizeMetrics, Integer> metricsFile = null;

    /**
     * Initialize the collector with input arguments;
     */
    @Override
    public void initialize(
            final InsertSizeMetricsArgumentCollection inputArgs,
            final SAMFileHeader samHeader,
            final List<Header> defaultHeaders)
    {
        this.inputArgs = inputArgs;
        collector.initialize(inputArgs, samHeader);
        metricsFile = new MetricsFile<>();
        if (defaultHeaders != null) {
            defaultHeaders.stream().forEach(h -> metricsFile.addHeader(h));
        }
    }

    /**
     * Return the read filters to be used for this collector.
     *
     * @return List of read filters for this collector
     */
    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return collector.getDefaultReadFilters();
    }

    @Override
    public void collectMetrics(
            final JavaRDD<GATKRead> filteredReads,
            final SAMFileHeader samHeader)
    {
        if (filteredReads.isEmpty()) {
            throw new GATKException("No valid reads found in input file.");
        }

        resultMetrics = filteredReads
                .mapPartitions(
                    it -> {
                        it.forEachRemaining(r -> collector.acceptRecord(r.convertToSAMRecord(samHeader), null));
                        return Arrays.asList(collector).iterator();
                    })
                .reduce(collector::combine);
    }

    @Override
    public void saveMetrics(final String inputName) {
        resultMetrics.finish(metricsFile, inputName);
    }

}
