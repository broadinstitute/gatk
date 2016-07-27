package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

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
        Assert.assertEquals(GATKProtectedMathUtils.interquartileRange(new double[]{1, 2, 3, 4, 5}), 3, 0.001);
        Assert.assertEquals(GATKProtectedMathUtils.interquartileRange(new double[]{4, 2, 1, 5, 3}), 3, 0.001);
        Assert.assertEquals(GATKProtectedMathUtils.interquartileRange(new double[]{1, 2, 3, 4, 5, 6, 7, 8, 9}), 5, 0.001);
    }

    @Test
    public void testMean() {
        Assert.assertEquals (GATKProtectedMathUtils.mean(0, 1, 2), 1, 1e-10);
        Assert.assertEquals (GATKProtectedMathUtils.mean(0, 1), 0.5, 1e-10);
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
        Assert.assertEquals(columnMeans, new double[]{2, 5, 8});
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

    @Test
    public void testLinspace() {
        Assert.assertEquals(GATKProtectedMathUtils.createEvenlySpacedPoints(-5, 1, 5), new double[] {-5.0000, -3.5000, -2.0000, -0.5000, 1.0000});
        Assert.assertEquals(GATKProtectedMathUtils.createEvenlySpacedPoints(-5, 1, 0), new double[] {});
        Assert.assertEquals(GATKProtectedMathUtils.createEvenlySpacedPoints(-5, 1, 1), new double[] {1.0});
        Assert.assertEquals(GATKProtectedMathUtils.createEvenlySpacedPoints(-5, 1, 2), new double[] {-5.0, 1.0});
        Assert.assertEquals(GATKProtectedMathUtils.createEvenlySpacedPoints(-5, -5, 2), new double[] {-5.0, -5.0});
        Assert.assertEquals(GATKProtectedMathUtils.createEvenlySpacedPoints(1e-10, 1, 2), new double[] {1e-10, 1});
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLinspaceWithNegativePoints() {
        GATKProtectedMathUtils.createEvenlySpacedPoints(-5, 1, -2);
    }



    @Test
    public void testRandomSelectFlatProbability() {
        final RandomGenerator rg = RandomGeneratorFactory.createRandomGenerator(new Random(13));
        final int NUM_SAMPLES = 1000;
        final List<Integer> choices = Arrays.asList(0,1,2);
        final List<Integer> result = IntStream.range(0, NUM_SAMPLES)
                .map(n -> GATKProtectedMathUtils.randomSelect(choices, j -> 1.0 / choices.size(), rg))
                .boxed()
                .collect(Collectors.toList());
        Assert.assertEquals(result.stream().filter(n -> n == 0).count(), NUM_SAMPLES / choices.size(), 50);
    }

    @Test
    public void testRandomSelect() {
        final RandomGenerator rg = RandomGeneratorFactory.createRandomGenerator(new Random(13));
        final int NUM_SAMPLES = 1000;
        final List<Integer> choices = Arrays.asList(-1,0,1);
        final List<Integer> result = IntStream.range(0, NUM_SAMPLES)
                .map(n -> GATKProtectedMathUtils.randomSelect(choices, j -> j*j/2.0, rg))
                .boxed()
                .collect(Collectors.toList());
        Assert.assertEquals(result.stream().filter(n -> n==0).count(), 0);
        Assert.assertEquals(result.stream().filter(n -> n == 1).count(), NUM_SAMPLES / 2, 50);
    }

    @Test
    public void testColumnSum() {
        final double[][] array = { {1,2,3}, {5,5,5}, {7,8,9}, {-15, 2, 12}};
        final double[] guess = GATKProtectedMathUtils.columnSums(new Array2DRowRealMatrix(array));
        final double[] gt = {-2, 17, 29};
        Assert.assertEquals(guess.length, 3);
        Assert.assertEquals(guess, gt);
    }

    @Test
    public void testColumnNaN() {
        final double[][] array = { {1,2, Double.NaN}, {5,5,5}, {7,8,9}, {-15, 2, 12}};
        final double[] guess = GATKProtectedMathUtils.columnSums(new Array2DRowRealMatrix(array));
        final double[] gt = {-2, 17, Double.NaN};
        Assert.assertEquals(guess.length, 3);
        Assert.assertEquals(guess[0], gt[0]);
        Assert.assertEquals(guess[1], gt[1]);
        Assert.assertTrue(Double.isNaN(guess[2]));
    }

    @Test
    public void testColumnInf() {
        final double[][] array = { {1,2, Double.POSITIVE_INFINITY}, {5,5,5}, {7,8,9}, {-15, 2, 12}};
        final double[] guess = GATKProtectedMathUtils.columnSums(new Array2DRowRealMatrix(array));
        final double[] gt = {-2, 17, Double.POSITIVE_INFINITY};
        Assert.assertEquals(guess.length, 3);
        Assert.assertEquals(guess[0], gt[0]);
        Assert.assertEquals(guess[1], gt[1]);
        Assert.assertTrue(Double.isInfinite(guess[2]));
    }

    @Test
    public void testColumnNegInf() {
        final double[][] array = { {1,2, Double.NEGATIVE_INFINITY}, {5,5,5}, {7,8,9}, {-15, 2, 12}};
        final double[] guess = GATKProtectedMathUtils.columnSums(new Array2DRowRealMatrix(array));
        final double[] gt = {-2, 17, Double.NEGATIVE_INFINITY};
        Assert.assertEquals(guess.length, 3);
        Assert.assertEquals(guess[0], gt[0]);
        Assert.assertEquals(guess[1], gt[1]);
        Assert.assertTrue(Double.isInfinite(guess[2]));
    }

    @Test
    public void testRowSum() {
        final double[][] array = { {1,2,3}, {5,5,5}, {7,8,9}, {-15, 2, 12}};
        final double[] guess = GATKProtectedMathUtils.rowSums(new Array2DRowRealMatrix(array));
        final double[] gt = {6, 15, 24, -1};
        Assert.assertEquals(guess.length, 4);
        Assert.assertEquals(guess, gt);
    }

    @Test
    public void testRowNaN() {
        final double[][] array = { {1,2, Double.NaN}, {5,5,5}, {7,8,9}, {-15, 2, 12}};
        final double[] guess = GATKProtectedMathUtils.rowSums(new Array2DRowRealMatrix(array));
        final double[] gt = {Double.NaN, 15, 24, -1};
        Assert.assertEquals(guess.length, 4);
        Assert.assertEquals(guess[1], gt[1]);
        Assert.assertEquals(guess[2], gt[2]);
        Assert.assertEquals(guess[3], gt[3]);
        Assert.assertTrue(Double.isNaN(guess[0]));
    }

    @Test
    public void testRowInf() {
        final double[][] array = { {1,2, Double.POSITIVE_INFINITY}, {5,5,5}, {7,8,9}, {-15, 2, 12}};
        final double[] guess = GATKProtectedMathUtils.rowSums(new Array2DRowRealMatrix(array));
        final double[] gt = {Double.POSITIVE_INFINITY, 15, 24, -1};
        Assert.assertEquals(guess.length, 4);
        Assert.assertEquals(guess[1], gt[1]);
        Assert.assertEquals(guess[2], gt[2]);
        Assert.assertEquals(guess[3], gt[3]);
        Assert.assertTrue(Double.isInfinite(guess[0]));
    }

    @Test
    public void testRowNegInf() {
        final double[][] array = { {1,2, Double.NEGATIVE_INFINITY}, {5,5,5}, {7,8,9}, {-15, 2, 12}};
        final double[] guess = GATKProtectedMathUtils.rowSums(new Array2DRowRealMatrix(array));
        final double[] gt = {Double.NEGATIVE_INFINITY, 15, 24, -1};
        Assert.assertEquals(guess.length, 4);
        Assert.assertEquals(guess[1], gt[1]);
        Assert.assertEquals(guess[2], gt[2]);
        Assert.assertEquals(guess[3], gt[3]);
        Assert.assertTrue(Double.isInfinite(guess[0]));
    }

    @Test
    public void testSum() {
        final double[] array = {-2, 17, 29};
        final double gt = 44;
        final double guess = MathUtils.sum(array);
        Assert.assertEquals(guess, gt);
    }

    @Test
    public void testSum3d() {
        final double[][][] array = {{{-2, 17}, {1, 1}}, {{1, 1}, {-2, 17}}, {{1, 1}, {1, 1}}, {{1, 1}, {1, 1}}};
        final double gt = 42;
        final double guess = GATKProtectedMathUtils.sum(array);
        Assert.assertEquals(guess, gt);
    }

    @Test
    public void testSum3dInf() {
        final double[][][] array = {{{Double.POSITIVE_INFINITY, 17}, {1, 1}}, {{1, 1}, {-2, 17}}, {{1, 1}, {1, 1}}, {{1, 1}, {1, 1}}};
        final double guess = GATKProtectedMathUtils.sum(array);
        Assert.assertTrue(Double.isInfinite(guess));
    }

    @Test
    public void testSum3dNaN() {
        final double[][][] array = {{{Double.NaN, 17}, {1, 1}}, {{1, 1}, {-2, 17}}, {{1, 1}, {1, 1}}, {{1, 1}, {1, 1}}};
        final double guess = GATKProtectedMathUtils.sum(array);
        Assert.assertTrue(Double.isNaN(guess));
    }

    @Test
    public void testSecondSmallestMinusSmallest() {
        Assert.assertEquals(GATKProtectedMathUtils.secondSmallestMinusSmallest(new int[0], -1), -1);
        Assert.assertEquals(GATKProtectedMathUtils.secondSmallestMinusSmallest(new int[] { 10 }, 3), 3);
        Assert.assertEquals(GATKProtectedMathUtils.secondSmallestMinusSmallest(new int[] { -10, -23, 3}, -1), 13);
        Assert.assertEquals(GATKProtectedMathUtils.secondSmallestMinusSmallest(new int[] { -10, -23, 3}, -1), 13);
    }

    @Test
    public void testMaxDifference() {
        final double[] array1 = {0.0, 1.0, 2.0};
        final double[] array2 = {-0.1, 1.05, 2.0};
        Assert.assertEquals(GATKProtectedMathUtils.maxDifference(array1, array2), 0.1, 1e-10);

        final double[] array3 = {0.0, 1.0, 2.0, 3.0};
        final double[] array4 = {0.0, 1.0, 2.0, 0.0};
        Assert.assertEquals(GATKProtectedMathUtils.maxDifference(array3, array4), 3.0, 1e-10);
    }
}