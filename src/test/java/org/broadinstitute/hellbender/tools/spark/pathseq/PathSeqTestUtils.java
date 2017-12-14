package org.broadinstitute.hellbender.tools.spark.pathseq;

public final class PathSeqTestUtils {

    /**
     * Returns true if a and b are within absoluteTolerance of each other.
     */
    public static boolean equalWithinTolerance(final double a, final double b, final double absoluteTolerance) {
        return Math.abs(a - b) < absoluteTolerance;
    }
}
