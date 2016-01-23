package org.broadinstitute.hellbender.tools.exome;

/**
 * Transformation performed to the evidence count.
 */
enum CoverageTransformation {

    /**
     * The final coverage figure is the absolute number of overlapping data.
     */
    NONE,

    /**
     * The final coverage figure is the average number of overlapping data per
     * base-pair in the target.
     */
    AVERAGE_PER_BP
}
