package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import htsjdk.samtools.metrics.MetricBase;

/**
 * High level metrics about the presence of  outward- and inward-facing pairs
 * within a SAM file generated with a jumping library, produced by
 * the CollectJumpingLibraryMetrics program and usually stored in a file with
 * the extension ".jump_metrics".
 */
public final class JumpingLibraryMetrics extends MetricBase {

    /**
     * The number of outward-facing pairs in the SAM file
     */
    public long JUMP_PAIRS;

    /**
     * The number of outward-facing pairs that are duplicates
     */
    public long JUMP_DUPLICATE_PAIRS;

    /**
     * The percentage of outward-facing pairs that are marked as duplicates
     */
    public double JUMP_DUPLICATE_PCT;

    /**
     * The estimated library size for outward-facing pairs
     */
    public long JUMP_LIBRARY_SIZE;

    /**
     * The mean insert size for outward-facing pairs
     */
    public double JUMP_MEAN_INSERT_SIZE;

    /**
     * The standard deviation on the insert size for outward-facing pairs
     */
    public double JUMP_STDEV_INSERT_SIZE;

    /**
     * The number of inward-facing pairs in the SAM file
     */
    public long NONJUMP_PAIRS;

    /**
     * The number of inward-facing pais that are duplicates
     */
    public long NONJUMP_DUPLICATE_PAIRS;

    /**
     * The percentage of inward-facing pairs that are marked as duplicates
     */
    public double NONJUMP_DUPLICATE_PCT;

    /**
     * The estimated library size for inward-facing pairs
     */
    public long NONJUMP_LIBRARY_SIZE;

    /**
     * The mean insert size for inward-facing pairs
     */
    public double NONJUMP_MEAN_INSERT_SIZE;

    /**
     * The standard deviation on the insert size for inward-facing pairs
     */
    public double NONJUMP_STDEV_INSERT_SIZE;

    /**
     * The number of pairs where either (a) the ends fall on different chromosomes or (b) the insert size
     * is greater than the maximum of 100000 or 2 times the mode of the insert size for outward-facing pairs.
     */
    public long CHIMERIC_PAIRS;

    /**
     * The number of fragments in the SAM file
     */
    public long FRAGMENTS;

    /**
     * The number of outward-facing pairs expressed as a percentage of the total of all outward facing pairs,
     * inward-facing pairs, and chimeric pairs.
     */
    public double PCT_JUMPS; 

    /**
     * The number of inward-facing pairs expressed as a percentage of the total of all outward facing pairs,
     * inward-facing pairs, and chimeric pairs.
     */
    public double PCT_NONJUMPS;

    /**
     * The number of chimeric pairs expressed as a percentage of the total of all outward facing pairs,
     * inward-facing pairs, and chimeric pairs.
     */
    public double PCT_CHIMERAS;
    
}
