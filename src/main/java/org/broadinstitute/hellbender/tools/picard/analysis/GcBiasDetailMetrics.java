package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.metrics.MetricBase;

/**
 * Class that holds detailed metrics about reads that fall within windows of a certain
 * GC bin on the reference genome.
 *
 * @author Tim Fennell
 */
public final class GcBiasDetailMetrics extends MetricBase {
    /** The G+C content of the reference sequence represented by this bin. Values are from 0% to 100% */
    public int GC;

    /** The number of windows on the reference genome that have this G+C content. */
    public int WINDOWS;

    /** The number of reads whose start position is at the start of a window of this GC. */
    public long READ_STARTS;

    /** The mean quality (determined via the error rate) of all bases of all reads that are assigned to windows of this GC. */
    public int MEAN_BASE_QUALITY;

    /**
     * The ration of "coverage" in this GC bin vs. the mean coverage of all GC bins. A number of
     * 1 represents mean coverage, a number less than one represents lower than mean coverage (e.g. 0.5
     * means half as much coverage as average) while a number greater than one represents higher than
     * mean coverage (e.g. 3.1 means this GC bin has 3.1 times more reads per window than average).
     */
    public double NORMALIZED_COVERAGE;

    /**
     * The radius of error bars in this bin based on the number of observations made. For example if
     * the normalized coverage is 0.75 and the error bar width is 0.1 then the error bars would be
     * drawn from 0.65 to 0.85. 
     */
    public double ERROR_BAR_WIDTH;

}
