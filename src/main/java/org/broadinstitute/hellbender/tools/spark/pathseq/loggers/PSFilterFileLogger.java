package org.broadinstitute.hellbender.tools.spark.pathseq.loggers;

import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.metrics.MetricsUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Logs filtering read counts to metrics file
 */
public final class PSFilterFileLogger implements PSFilterLogger {

    private final PSFilterMetrics metrics;
    private final MetricsFile<PSFilterMetrics, Long> metricsFile;
    private final String metricsOutputPath;

    public PSFilterFileLogger(final MetricsFile<PSFilterMetrics, Long> metricsFile, final String metricsOutputPath) {
        Utils.nonNull(metricsFile, "Filter logger parameter metricsFile cannot be null");
        Utils.nonNull(metricsOutputPath, "Filter logger parameter metricsOutputPath cannot be null");
        this.metrics = new PSFilterMetrics();
        this.metricsFile = metricsFile;
        this.metricsOutputPath = metricsOutputPath;
    }

    @Override
    public void logPrimaryReads(final JavaRDD<GATKRead> reads) {
        Utils.nonNull(reads, "Filter logging parameter reads cannot be null");
        metrics.PRIMARY_READS = reads.count();
    }

    @Override
    public void logReadsAfterPrealignedHostFilter(final JavaRDD<GATKRead> reads) {
        Utils.nonNull(reads, "Filter logging parameter reads cannot be null");
        metrics.READS_AFTER_PREALIGNED_HOST_FILTER = reads.count();
    }

    @Override
    public void logReadsAfterQualityFilter(final JavaRDD<GATKRead> reads) {
        Utils.nonNull(reads, "Filter logging parameter reads cannot be null");
        metrics.READS_AFTER_QUALITY_AND_COMPLEXITY_FILTER = reads.count();
    }

    @Override
    public void logReadsAfterHostFilter(final JavaRDD<GATKRead> reads) {
        Utils.nonNull(reads, "Filter logging parameter reads cannot be null");
        metrics.READS_AFTER_HOST_FILTER = reads.count();
    }

    @Override
    public void logReadsAfterDeduplication(final JavaRDD<GATKRead> reads) {
        Utils.nonNull(reads, "Filter logging parameter reads cannot be null");
        metrics.READS_AFTER_DEDUPLICATION = reads.count();
    }

    @Override
    public void logFinalPairedReads(final JavaRDD<GATKRead> reads) {
        Utils.nonNull(reads, "Filter logging parameter reads cannot be null");
        metrics.FINAL_PAIRED_READS = reads.count();
    }

    @Override
    public void close() {
        metrics.computeDerivedMetrics();
        metricsFile.addMetric(metrics);
        MetricsUtils.saveMetrics(metricsFile, metricsOutputPath);
    }
}
