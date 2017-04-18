package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;
import org.testng.Assert;

/**
 * This class provides useful assertions about approximate equality of various mathematical objects
 * (e.g. tensors, matrices, vectors, etc) and/or conditions imposed on them (e.g. symmetry, positive-definiteness, etc).
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class MathObjectAsserts {
    public static final double DEFAULT_RELATIVE_TOLERANCE = 1e-10;
    public static final double DEFAULT_ABSOLUTE_TOLERANCE = 1e-10;

    private MathObjectAsserts() {}

    public static void assertRealMatrixEquals(final RealMatrix actual, final RealMatrix expected,
                                              final double relativeTolerance, final double absoluteTolerance) {
        Assert.assertNotNull(expected);
        Assert.assertNotNull(actual);
        Assert.assertEquals(expected.getRowDimension(), actual.getRowDimension());
        Assert.assertEquals(expected.getColumnDimension(), actual.getColumnDimension());
        final int rowDimension = expected.getRowDimension();
        final int colDimension = expected.getColumnDimension();
        for (int i = 0; i < rowDimension; i++) {
            for (int j = 0; j < colDimension; j++) {
                assertDoubleEquals(expected.getEntry(i, j), actual.getEntry(i, j), relativeTolerance, absoluteTolerance,
                        String.format("Different entries at (%d, %d)", i, j));
            }
        }
        Assert.assertEquals(1.0, 1.0, 1.0);
    }

    public static void assertRealMatrixEquals(final RealMatrix actual, final RealMatrix expected,
                                              final double relativeTolerance) {
        assertRealMatrixEquals(actual, expected, relativeTolerance, DEFAULT_ABSOLUTE_TOLERANCE);
    }

    public static void assertRealMatrixEquals(final RealMatrix actual, final RealMatrix expected) {
        assertRealMatrixEquals(actual, expected, DEFAULT_RELATIVE_TOLERANCE, DEFAULT_ABSOLUTE_TOLERANCE);
    }

    public static void assertDoubleEquals(final double actual, final double expected,
                                          final double relativeTolerance, final double absoluteTolerance,
                                          final String message) {
        if (FastMath.abs(expected - actual) < absoluteTolerance)
            return;
        final double relativeError;
        if (FastMath.abs(actual) > FastMath.abs(expected)) {
            relativeError = FastMath.abs((expected - actual) / actual);
        } else {
            relativeError = FastMath.abs((expected - actual) / expected);
        }
        if (relativeError <= relativeTolerance)
            return;
        throw new AssertionError(String.format("%s; expected: %f, actual: %f", message, expected, actual));
    }

    public static void assertDoubleEquals(final double actual, final double expected,
                                          final double relativeTolerance, final double absoluteTolerance) {
        assertDoubleEquals(expected, actual, relativeTolerance, absoluteTolerance, "Unequal floating point numbers");
    }

    public static void assertDoubleEquals(final double actual, final double expected,
                                          final double relativeTolerance) {
        assertDoubleEquals(expected, actual, relativeTolerance, DEFAULT_ABSOLUTE_TOLERANCE);
    }

    public static void assertDoubleEquals(final double actual, final double expected) {
        assertDoubleEquals(expected, actual, DEFAULT_RELATIVE_TOLERANCE, DEFAULT_ABSOLUTE_TOLERANCE);
    }
}
