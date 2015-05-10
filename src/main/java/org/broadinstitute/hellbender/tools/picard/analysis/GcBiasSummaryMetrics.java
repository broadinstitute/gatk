package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.metrics.MetricBase;

/**
 * High level metrics that capture how biased the coverage in a certain lane is.
 *
 * @author Tim Fennell
 */
public final class GcBiasSummaryMetrics extends MetricBase {
    /** The window size on the genome used to calculate the GC of the sequence. */
    public int WINDOW_SIZE;

    /** The total number of clusters that were seen in the gc bias calculation. */
    public int TOTAL_CLUSTERS;

    /** The total number of aligned reads used to compute the gc bias metrics. */
    public int ALIGNED_READS;

    /**
     * Illumina-style AT dropout metric.  Calculated by taking each GC bin independently and calculating
     * (%ref_at_gc - %reads_at_gc) and summing all positive values for GC=[0..50].
     */
    public double AT_DROPOUT;

    /**
     * Illumina-style GC dropout metric.  Calculated by taking each GC bin independently and calculating
     * (%ref_at_gc - %reads_at_gc) and summing all positive values for GC=[50..100].
     */
    public double GC_DROPOUT;

}
