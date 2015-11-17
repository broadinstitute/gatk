package org.broadinstitute.hellbender.tools.exome;

/**
 * The piece of evidence used to quantify coverage.
 */
enum CoverageDatum {

    /**
     * Coverage is measured by the number of overlapping reads.
     */
    READ,

    /**
     * Coverage is measured by the number of overlapping base calls.
     */
    BASE_CALL,
}
