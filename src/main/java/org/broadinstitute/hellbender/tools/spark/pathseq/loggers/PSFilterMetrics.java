package org.broadinstitute.hellbender.tools.spark.pathseq.loggers;

import htsjdk.samtools.metrics.MetricBase;

/**
 * Metrics that are calculated during the PathSeq filter
 */
public final class PSFilterMetrics extends MetricBase {

    /**
     * Number of input reads that are not secondary or supplementary alignments
     */
    public Long PRIMARY_READS;

    /**
     * The number of reads after filtering prealigned reads
     */
    public Long READS_AFTER_PREALIGNED_HOST_FILTER;

    /**
     * The number of reads after low-quality and low-complexity filtering
     */
    public Long READS_AFTER_QUALITY_AND_COMPLEXITY_FILTER;

    /**
     * The number of reads after host read filtering
     */
    public Long READS_AFTER_HOST_FILTER;

    /**
     * The number of reads after deduplication
     */
    public Long READS_AFTER_DEDUPLICATION;

    /**
     * The number of non-host reads that have mates
     */
    public Long FINAL_PAIRED_READS;

    /**
     * The number of non-host unpaired reads
     */
    public Long FINAL_UNPAIRED_READS;

    /**
     * The number of non-host reads in total
     */
    public Long FINAL_TOTAL_READS;

    /**
     * The number of reads filtered because they were low-quality or low-complexity
     */
    public Long LOW_QUALITY_OR_LOW_COMPLEXITY_READS_FILTERED;

    /**
     * The number of subtracted reads identified as host
     */
    public Long HOST_READS_FILTERED;

    /**
     * The number of filtered duplicate reads
     */
    public Long DUPLICATE_READS_FILTERED;

    /**
     * Computes number of reads filtered at each step and total final reads.
     */
    public void computeDerivedMetrics() {
        if (PRIMARY_READS == null || READS_AFTER_PREALIGNED_HOST_FILTER == null
                || READS_AFTER_QUALITY_AND_COMPLEXITY_FILTER == null || READS_AFTER_HOST_FILTER == null
                || READS_AFTER_DEDUPLICATION == null || FINAL_PAIRED_READS == null) {
            throw new IllegalStateException("Cannot compute metrics if primary, pre-aligned host, quality, host, duplicate, or final paired read counts are not initialized");
        }
        FINAL_TOTAL_READS = READS_AFTER_DEDUPLICATION;
        FINAL_UNPAIRED_READS = FINAL_TOTAL_READS - FINAL_PAIRED_READS;
        LOW_QUALITY_OR_LOW_COMPLEXITY_READS_FILTERED = READS_AFTER_PREALIGNED_HOST_FILTER - READS_AFTER_QUALITY_AND_COMPLEXITY_FILTER;
        HOST_READS_FILTERED = PRIMARY_READS - READS_AFTER_PREALIGNED_HOST_FILTER + READS_AFTER_QUALITY_AND_COMPLEXITY_FILTER - READS_AFTER_HOST_FILTER;
        DUPLICATE_READS_FILTERED = READS_AFTER_HOST_FILTER - READS_AFTER_DEDUPLICATION;
    }
}
