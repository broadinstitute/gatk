package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.tools.pon.PoNTestUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class MatrixSummaryUtilsUnitTest extends BaseTest {

    @Test(dataProvider = "rowAndColumnMedianTestData")
    public void testGetRowAndColumnMedians(double[][] testdata, double[] gtRowMedian, double[] gtColumnMedian) {
        final RealMatrix tmp = new Array2DRowRealMatrix(testdata);
        PoNTestUtils.assertEqualsDoubleArrays(MatrixSummaryUtils.getRowMedians(tmp), gtRowMedian);
        PoNTestUtils.assertEqualsDoubleArrays(MatrixSummaryUtils.getColumnMedians(tmp), gtColumnMedian);
    }

    @Test(dataProvider = "rowVarianceTestData")
    public void testGetRowVariances(double[][] testdata, double[] gtRowVariances) {
        final RealMatrix tmp = new Array2DRowRealMatrix(testdata);
        PoNTestUtils.assertEqualsDoubleArrays(MatrixSummaryUtils.getRowVariances(tmp), gtRowVariances);
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
}
