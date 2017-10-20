package org.broadinstitute.hellbender.tools.spark.pathseq.loggers;

import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.metrics.MetricsUtils;
import org.broadinstitute.hellbender.tools.spark.pathseq.PSScorer;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Logs number of mapped and unmapped reads to metrics file
 */
public final class PSScoreFileLogger implements PSScoreLogger {

    private final PSScoreMetrics metrics;
    private final MetricsFile<PSScoreMetrics, Long> metricsFile;
    private final String metricsOutputPath;

    public PSScoreFileLogger(final MetricsFile<PSScoreMetrics, Long> metricsFile, final String metricsOutputPath) {
        this.metrics = new PSScoreMetrics();
        this.metricsFile = metricsFile;
        this.metricsOutputPath = metricsOutputPath;
    }

    public void logReadCounts(final JavaRDD<GATKRead> reads) {
        final long numReads = reads.count();
        final long numMappedReads = reads.filter(read -> read.hasAttribute(PSScorer.HITS_TAG)).count();
        final long numUnmappedReads = numReads - numMappedReads;
        metrics.UNMAPPED_READS = numUnmappedReads;
        metrics.MAPPED_READS = numMappedReads;
    }

    public void close() {
        metricsFile.addMetric(metrics);
        MetricsUtils.saveMetrics(metricsFile, metricsOutputPath);
    }
}
