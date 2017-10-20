package org.broadinstitute.hellbender.tools.spark.pathseq.loggers;

import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.metrics.MetricsUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Logs filtering read counts to metrics file
 */
public final class PSFilterFileLogger implements PSFilterLogger {

    private final PSFilterMetrics metrics;
    private final MetricsFile<PSFilterMetrics, Long> metricsFile;
    private final String metricsOutputPath;

    public PSFilterFileLogger(final MetricsFile<PSFilterMetrics, Long> metricsFile, final String metricsOutputPath) {
        this.metrics = new PSFilterMetrics();
        this.metricsFile = metricsFile;
        this.metricsOutputPath = metricsOutputPath;
    }

    public void logPrimaryReads(final JavaRDD<GATKRead> reads) { metrics.PRIMARY_READS = reads.count(); }
    public void logReadsAfterPrealignedHostFilter(final JavaRDD<GATKRead> reads) { metrics.READS_AFTER_PREALIGNED_HOST_FILTER = reads.count(); }
    public void logReadsAfterQualityFilter(final JavaRDD<GATKRead> reads) { metrics.READS_AFTER_QUALITY_AND_COMPLEXITY_FILTER = reads.count(); }
    public void logReadsAfterHostFilter(final JavaRDD<GATKRead> reads) { metrics.READS_AFTER_HOST_FILTER = reads.count(); }
    public void logReadsAfterDeduplication(final JavaRDD<GATKRead> reads) { metrics.READS_AFTER_DEDUPLICATION = reads.count(); }
    public void logFinalPairedReads(final JavaRDD<GATKRead> reads) { metrics.FINAL_PAIRED_READS = reads.count(); }

    public void close() {
        metrics.computeDerivedMetrics();
        metricsFile.addMetric(metrics);
        MetricsUtils.saveMetrics(metricsFile, metricsOutputPath);
    }
}
