package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class MatrixSummaryUtilsUnitTest extends GATKBaseTest {
    private static final double DOUBLE_ARRAY_TOLERANCE = 1E-8;

    @Test(dataProvider = "rowAndColumnMedianTestData")
    public void testGetRowAndColumnMedians(double[][] testdata, double[] gtRowMedian, double[] gtColumnMedian) {
        final RealMatrix tmp = new Array2DRowRealMatrix(testdata);
        assertEqualsDoubleArrays(MatrixSummaryUtils.getRowMedians(tmp), gtRowMedian, DOUBLE_ARRAY_TOLERANCE);
        assertEqualsDoubleArrays(MatrixSummaryUtils.getColumnMedians(tmp), gtColumnMedian, DOUBLE_ARRAY_TOLERANCE);
    }

    @Test(dataProvider = "rowVarianceTestData")
    public void testGetRowVariances(double[][] testdata, double[] gtRowVariances) {
        final RealMatrix tmp = new Array2DRowRealMatrix(testdata);
        assertEqualsDoubleArrays(MatrixSummaryUtils.getRowVariances(tmp), gtRowVariances, DOUBLE_ARRAY_TOLERANCE);
    }

    @DataProvider(name = "rowAndColumnMedianTestData")
    public Object[][] createRowAndColumnMedianTestData() {
        return new Object[][] {
                {
                        new double[][]{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}},
                        new double[]{2.0, 5.0, 8.0},
                        new double[]{4.0, 5.0, 6.0}
                },
                {
                        new double[][]{{Double.NaN, Double.NaN, 4.0}, {1.0, 2.0, 3.0}},
                        new double[] {4.0, 2.0},
                        new double[] {1.0, 2.0, 3.5},
                }
        };
    }

    @DataProvider(name = "rowVarianceTestData")
    public Object[][] createRowVarianceTestData() {
        return new Object[][] {
                {
                        new double[][]{{1, 2, 3}, {4, 5, 6,}, {7, 8, 9.1}, {10, 11.25, 15}},
                        new double[]{1.0, 1.0, 1.10333333333333, 6.77083333333333},
                },
                {
                        new double[][]{{Double.NaN, 3.0, 4.0, 5.0}, {10.0, 2.0, 3.0, 4.0}},
                        new double[] {Double.NaN, 12.9166666666667},
                },
                {
                        new double[][]{{Double.NaN, Double.NaN, 4.0}, {1.0, 2.0, 3.0}},
                        new double[] {Double.NaN, 1.0},
                }
        };
    }

    /**
     * Test whether two double arrays are equal.  For this method NaN is considered to equal NaN
     *
     * @param actual never {@code null}
     * @param gt never {@code null}
     */
    private static void assertEqualsDoubleArrays(final double[] actual, final double[] gt, final double tolerance) {
        Assert.assertEquals(actual.length, gt.length);
        for (int i = 0; i < actual.length; i++) {
            if (Double.isNaN(gt[i])) {
                Assert.assertTrue(Double.isNaN(actual[i]));
            } else {
                Assert.assertEquals(actual[i], gt[i], tolerance, String.format("Arrays were not equal (within tolerance %s) at index %d.", Double.toString(tolerance), i));
            }
        }
    }
}
