package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * Created by davidben on 1/22/16.
 */
public class GATKProtectedMathUtilsTest {
    @Test
    public void naturalLogSumNaturalLogTest() {
        final double[] realSpaceValues = {0.1, 0.2, 1.4, 5.9};
        final double[] logSpaceValues= Arrays.stream(realSpaceValues).map(Math::log).toArray();

        final double expected = Math.log(Arrays.stream(realSpaceValues).sum());
        final double actual = GATKProtectedMathUtils.logSumExp(logSpaceValues);
        Assert.assertEquals(actual, expected, 1e-10);
    }

	@Test
    public void meanTest() {
        Assert.assertEquals(GATKProtectedMathUtils.mean(1,2,3), 2, 1e-10);
        Assert.assertEquals(GATKProtectedMathUtils.mean(1,2), 1.5, 1e-10);
        Assert.assertEquals(GATKProtectedMathUtils.mean(new double[] {1,2}), 1.5, 1e-10);
    }

    @Test
    public void interquartileRangeTest() {
        Assert.assertEquals(GATKProtectedMathUtils.interquartileRange(new double[] {1,2,3,4,5}), 3, 0.001);
        Assert.assertEquals(GATKProtectedMathUtils.interquartileRange(new double[] {4,2,1,5,3}), 3, 0.001);
        Assert.assertEquals(GATKProtectedMathUtils.interquartileRange(new double[] {1,2,3,4,5,6,7,8,9}), 5, 0.001);
    }

    @Test
    public void testMean() {
        Assert.assertEquals (GATKProtectedMathUtils.mean(0,1,2), 1, 1e-10);
        Assert.assertEquals (GATKProtectedMathUtils.mean(0,1), 0.5, 1e-10);
    }

    @Test
    public void testRowMeans() {
        final double[][] array = { {1,2,3}, {4,5,6}, {7,8,9}};
        final double[] rowMeans = GATKProtectedMathUtils.rowMeans(new Array2DRowRealMatrix(array));
        Assert.assertEquals(rowMeans, new double[] {2,5,8});
    }

    @Test
    public void testRowVariances() {
        final double[][] array = { {1,2,3}, {5,5,5}, {7,8,9}};
        final double[] rowVariances = GATKProtectedMathUtils.rowVariances(new Array2DRowRealMatrix(array));
        Assert.assertEquals(rowVariances[0], 1, 1e-8);
        Assert.assertEquals(rowVariances[1], 0, 1e-8);
        Assert.assertEquals(rowVariances[2], 1, 1e-8);
    }

    @Test
    public void testRowStdDevs() {
        final double[][] array = { {1,2,3}, {5,5,5}, {7,8,9}, {-15, 2, 12}};
        final double[] rowStdDevs = GATKProtectedMathUtils.rowStdDevs(new Array2DRowRealMatrix(array));
        Assert.assertEquals(rowStdDevs[0], 1, 1e-8);
        Assert.assertEquals(rowStdDevs[1], 0, 1e-8);
        Assert.assertEquals(rowStdDevs[2], 1, 1e-8);
        Assert.assertEquals(rowStdDevs[3], 13.65039682, 1e-8);
    }

    @Test
    public void testStdDevsOnAList() {
        Assert.assertEquals(GATKProtectedMathUtils.stdDev(Arrays.asList(1, 2, 3)), 1, 1e-8);
        Assert.assertEquals(GATKProtectedMathUtils.stdDev(Arrays.asList(5, 5, 5)), 0, 1e-8);
        Assert.assertEquals(GATKProtectedMathUtils.stdDev(Arrays.asList(7, 8, 9)), 1, 1e-8);
        Assert.assertEquals(GATKProtectedMathUtils.stdDev(Arrays.asList(-15, 2, 12)), 13.65039682, 1e-8);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testStdDevsOnAListWithNull() {
        GATKProtectedMathUtils.stdDev(Arrays.asList(1, null, 3));
    }

    @Test
    public void testStdDevsOnAnArray() {
        Assert.assertEquals(GATKProtectedMathUtils.stdDev(1, 2, 3), 1, 1e-8);
        Assert.assertEquals(GATKProtectedMathUtils.stdDev(5, 5, 5), 0, 1e-8);
        Assert.assertEquals(GATKProtectedMathUtils.stdDev(7, 8, 9), 1, 1e-8);
        Assert.assertEquals(GATKProtectedMathUtils.stdDev(-15, 2, 12), 13.65039682, 1e-8);
    }

    @Test
    public void testColumnMeans() {
        final double[][] array = { {1,2,3}, {4,5,6}, {7,8,9}};
        final double[] columnMeans = GATKProtectedMathUtils.columnMeans(new Array2DRowRealMatrix(array).transpose());
        Assert.assertEquals(columnMeans, new double[] {2,5,8});
    }

    @Test
    public void testColumnVariances() {
        final double[][] array = { {1,2,3}, {5,5,5}, {7,8,9}};
        final double[] columnVariances = GATKProtectedMathUtils.columnVariances(new Array2DRowRealMatrix(array).transpose());
        Assert.assertEquals(columnVariances[0], 1, 1e-8);
        Assert.assertEquals(columnVariances[1], 0, 1e-8);
        Assert.assertEquals(columnVariances[2], 1, 1e-8);
    }

    @Test
    public void testColumnStdDevs() {
        final double[][] array = { {1,2,3}, {5,5,5}, {7,8,9}, {-15, 2, 12}};
        final double[] columnStdDevs = GATKProtectedMathUtils.columnStdDevs(new Array2DRowRealMatrix(array).transpose());
        Assert.assertEquals(columnStdDevs[0], 1, 1e-8);
        Assert.assertEquals(columnStdDevs[1], 0, 1e-8);
        Assert.assertEquals(columnStdDevs[2], 1, 1e-8);
        Assert.assertEquals(columnStdDevs[3], 13.65039682, 1e-8);
    }

}