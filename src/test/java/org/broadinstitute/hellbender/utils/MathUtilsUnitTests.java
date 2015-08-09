package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Random;

import static java.lang.Math.log10;

/**
 * Basic unit test for MathUtils
 */
public final class MathUtilsUnitTests extends BaseTest {

    private final static Logger logger = LogManager.getLogger(MathUtilsUnitTests.class);

    @BeforeClass
    public void init() {
    }

    @Test
    public void testGoodProbability(){
        for (final double good : Arrays.asList(0.1, 1.0)) {
            Assert.assertTrue(MathUtils.goodProbability(good));
            Assert.assertTrue(MathUtils.goodLogProbability(Math.log(good)));
            Assert.assertTrue(MathUtils.goodLogProbability(Math.log(good), true));
            Assert.assertTrue(MathUtils.goodLogProbability(Math.log(good), false));
        }

        Assert.assertTrue(MathUtils.goodProbability(0.0));
        Assert.assertTrue(MathUtils.goodLogProbability(Double.NEGATIVE_INFINITY));
        Assert.assertTrue(MathUtils.goodLogProbability(Double.NEGATIVE_INFINITY, true));
        Assert.assertFalse(MathUtils.goodLogProbability(Double.NEGATIVE_INFINITY, false));

        for (final double bad : Arrays.asList(-1.0, 2.0, Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY)) {
            Assert.assertFalse(MathUtils.goodProbability(bad));
        }
    }

    @Test
    public void testLogOneMinusX(){
        Assert.assertEquals(-0.10536051565, MathUtils.logOneMinusX(0.1), 1e-6); //result from Wolfram Alpha
        Assert.assertEquals(0.0, MathUtils.logOneMinusX(0.0), 1e-6);
        Assert.assertEquals(Double.NEGATIVE_INFINITY, MathUtils.logOneMinusX(1.0), 1e-6);
    }

    @Test
    public void testBinomialProbability() {
        // results from Wolfram Alpha
        Assert.assertEquals(MathUtils.binomialProbability(3, 2, 0.5), 0.375, 0.0001);
        Assert.assertEquals(MathUtils.binomialProbability(100, 10, 0.5), 1.365543e-17, 1e-18);
        Assert.assertEquals(MathUtils.binomialProbability(217, 73, 0.02), 4.521904e-67, 1e-68);
        Assert.assertEquals(MathUtils.binomialProbability(300, 100, 0.02), 9.27097e-91, 1e-92);
        Assert.assertEquals(MathUtils.binomialProbability(300, 150, 0.98), 6.462892e-168, 1e-169);
        Assert.assertEquals(MathUtils.binomialProbability(300, 120, 0.98), 3.090054e-221, 1e-222);
        Assert.assertEquals(MathUtils.binomialProbability(300, 112, 0.98), 2.34763e-236, 1e-237);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBinomialProbabilityError() {
        Assert.assertEquals(MathUtils.binomialProbability(3, 2, 1.5), 0.375, 0.0001);
    }

    @Test
    public void testLogBinomialCoefficient() {
        // note that we can test the binomial coefficient calculation indirectly via Newton's identity
        // (1+z)^m = sum (m choose k)z^k
        double[] z_vals = new double[]{0.999, 0.9, 0.8, 0.5, 0.2, 0.01, 0.0001};
        int[] exponent = new int[]{5, 15, 25, 50, 100};
        for (double z : z_vals) {
            double logz = Math.log(z);
            for (int exp : exponent) {
                double expected_log = exp * Math.log(1 + z);
                double[] newtonArray_log = new double[1 + exp];
                for (int k = 0; k <= exp; k++) {
                    newtonArray_log[k] = MathUtils.logBinomialCoefficient(exp, k) + k * logz;
                }
                Assert.assertEquals(MathUtils.logSumLog(newtonArray_log), expected_log, 1e-6);
            }
        }

        // results from Wolfram Alpha
        Assert.assertEquals(MathUtils.logBinomialCoefficient(4, 2), 1.7917595, 1e-6);
        Assert.assertEquals(MathUtils.logBinomialCoefficient(10, 3), 4.7874917, 1e-6);
        Assert.assertEquals(MathUtils.logBinomialCoefficient(103928, 119), 921.5305037, 1e-4);
    }


    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogBinomialCoefficientErrorN() {
        Assert.assertEquals(MathUtils.logBinomialCoefficient(-1, 1), 0.0, 1e-6);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogBinomialCoefficientErrorK() {
        Assert.assertEquals(MathUtils.logBinomialCoefficient(1, -1), 0.0, 1e-6);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLogBinomialCoefficientErrorKmoreThanN() {
        Assert.assertEquals(MathUtils.logBinomialCoefficient(1, 2), 0.0, 1e-6);
    }

    @Test
    public void testApproximateLogSumLog() {
        final double requiredPrecision = 1E-4;

        Assert.assertEquals(MathUtils.approximateLogSumLog(new double[]{0.0, 0.0, 0.0}), Math.log(3), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLogSumLog(new double[]{0.0, 0.0, 0.0}, 0), 0.0, requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLogSumLog(new double[]{0.0, 0.0, 0.0}, 3), Math.log(3), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLogSumLog(new double[]{0.0, 0.0, 0.0}, 2), Math.log(2), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLogSumLog(new double[]{0.0, 0.0, 0.0}, 1), 0.0, requiredPrecision);

        Assert.assertEquals(MathUtils.approximateLogSumLog(new double[]{Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY}), Double.NEGATIVE_INFINITY);

        final Random random = new Random(13);
        for (int j = 0; j < 5; j++) {
            for (int i = 0; i < 5; i++) {
                final double a = (1 + 3 * j) * random.nextGaussian();
                final double b = (1 + 3 * j) * random.nextGaussian();
                final double c = (1 + 3 * j) * random.nextGaussian();

                Assert.assertEquals(MathUtils.approximateLogSumLog(new double[]{a}), a, requiredPrecision);
                Assert.assertEquals(MathUtils.approximateLogSumLog(new double[]{a, Double.NEGATIVE_INFINITY}), a, requiredPrecision);
                Assert.assertEquals(MathUtils.approximateLogSumLog(new double[]{a, b}), Math.log(Math.exp(a) + Math.exp(b)), requiredPrecision);
                Assert.assertEquals(MathUtils.approximateLogSumLog(a, b), Math.log(Math.exp(a) + Math.exp(b)), requiredPrecision);
                Assert.assertEquals(MathUtils.approximateLogSumLog(b, a), Math.log(Math.exp(a) + Math.exp(b)), requiredPrecision);
                Assert.assertEquals(MathUtils.approximateLogSumLog(new double[]{a, b, c}), Math.log(Math.exp(a) + Math.exp(b) + Math.exp(c)), requiredPrecision);
                Assert.assertEquals(MathUtils.approximateLogSumLog(a, b, c), Math.log(Math.exp(a) + Math.exp(b) + Math.exp(c)), requiredPrecision);
            }
        }
    }

    @Test
    public void testLogSumLog() {
        final double requiredPrecision = 1E-14;

        Assert.assertEquals(MathUtils.logSumLog(new double[]{0.0, 0.0, 0.0}), Math.log(3), requiredPrecision);
        Assert.assertEquals(MathUtils.logSumLog(new double[]{0.0, 0.0, 0.0}, 0), Math.log(3), requiredPrecision);
        Assert.assertEquals(MathUtils.logSumLog(new double[]{0.0, 0.0, 0.0}, 0, 3), Math.log(3), requiredPrecision);
        Assert.assertEquals(MathUtils.logSumLog(new double[]{0.0, 0.0, 0.0}, 0, 2), Math.log(2), requiredPrecision);
        Assert.assertEquals(MathUtils.logSumLog(new double[]{0.0, 0.0, 0.0}, 0, 1), 0.0, requiredPrecision);

        final Random random = new Random(13);
        for (int j = 0; j < 5; j++) {
            for (int i = 0; i < 5; i++) {
                final double a = (1 + 3 * j) * random.nextGaussian();
                final double b = (1 + 3 * j) * random.nextGaussian();
                final double c = (1 + 3 * j) * random.nextGaussian();

                Assert.assertEquals(MathUtils.logSumLog(new double[]{a}), a, requiredPrecision);
                Assert.assertEquals(MathUtils.logSumLog(new double[]{a, Double.NEGATIVE_INFINITY}), a, requiredPrecision);
                Assert.assertEquals(MathUtils.logSumLog(new double[]{a, b}), Math.log(Math.exp(a) + Math.exp(b)), requiredPrecision);
                Assert.assertEquals(MathUtils.logSumLog(new double[]{a, b, c}), Math.log(Math.exp(a) + Math.exp(b) + Math.exp(c)), requiredPrecision);
            }
        }
    }

    @Test
    public void testlogsumLogEdgeCases(){
        final double requiredPrecision = 1E-14;

        Assert.assertEquals(MathUtils.logSumLog(new double[]{3.0, 2.0, 1.0}, 2, 1), Double.NEGATIVE_INFINITY);
        Assert.assertEquals(MathUtils.logSumLog(new double[]{Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY}), Double.NEGATIVE_INFINITY);
        Assert.assertEquals(MathUtils.logSumLog(new double[]{3.0, 2.0, 1.0}), Math.log(Math.exp(3.0) + Math.exp(2.0) + Math.exp(1.0)), requiredPrecision);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testlogsumLogErrorNaN(){
        MathUtils.logSumLog(new double[]{Double.NaN, 1.0});
    }

    @Test
    public void testNormalize(){
        final double error = 1e-6;

        final double[] normalized = MathUtils.normalizeFromLog(new double[]{Math.log(3.0), Math.log(2.0), Math.log(1.0)});
        final double[] normalizedLog = MathUtils.normalizeFromLog(new double[]{Math.log(3.0), Math.log(2.0), Math.log(1.0)}, true);
        final double[] normalizedLogInLog = MathUtils.normalizeFromLog(new double[]{3.0, 2.0, 1.0}, true, true);
        final double[] normalizedExpected    = {3.0/6.0, 2.0/6.0, 1.0/6.0};  //sum of those is 1
        final double[] normalizedLogExpected = {Math.log(3.0/6.0), Math.log(2.0/6.0), Math.log(1.0/6.0)}; //sum of exp of those is 1
        final double[] normalizedLogInLogExpected = {0.0, -1.0, -2.0};

        assertEqualsDoubleArray(normalizedExpected, normalized, error);
        assertEqualsDoubleArray(normalizedLogExpected, normalizedLog, error);
        assertEqualsDoubleArray(normalizedLogInLogExpected, normalizedLogInLog, error);
    }

    @Test
    public void testLogFactorial() {
        logger.warn("Executing testLogFactorial");

        // results from Wolfram Alpha
        Assert.assertEquals(MathUtils.logFactorial(4), 3.1780538, 1e-6);
        Assert.assertEquals(MathUtils.logFactorial(10), 15.1044125, 1e-6);
        Assert.assertEquals(MathUtils.logFactorial(200), 863.2319872, 1e-3);
        Assert.assertEquals(MathUtils.logFactorial(12342), 103934.6907023, 1e-1);

        double logFactorial_small = 0;
        double logFactorial_middle = 863.2319872;
        double logFactorial_large = 103934.6907023;
        int small_start = 1;
        int med_start = 200;
        int large_start = 12342;
        for ( int i = 1; i < 1000; i++ ) {
            logFactorial_small += Math.log(i + small_start);
            logFactorial_middle += Math.log(i + med_start);
            logFactorial_large += Math.log(i + large_start);

            Assert.assertEquals(MathUtils.logFactorial(small_start+i),logFactorial_small,1e-6);
            Assert.assertEquals(MathUtils.logFactorial(med_start+i),logFactorial_middle,1e-3);
            Assert.assertEquals(MathUtils.logFactorial(large_start+i),logFactorial_large,1e-1);
        }
    }

    @Test
    public void testCovarianceDivergences() {
        logger.warn("Executing testCovarianceDivergences");
        //two symmetric positive-definite matrices
        double[][] cov1 = { {5, 2, 3},
                            {2, 7, 5},
                            {3, 5, 6}};

        double[][] cov2 = { {11, 3, 3},
                            {3, 7, 5},
                            {3, 5, 13}};

        RealMatrix mat1 = new Array2DRowRealMatrix(cov1);
        RealMatrix mat2 = new Array2DRowRealMatrix(cov2);

        Assert.assertEquals(MathUtils.covarianceKLDivergence(mat1, mat2), 3.65393, 1e-4);   //from Mathematica
        Assert.assertEquals(MathUtils.covarianceGeodesicDistance(mat1, mat2), 1.86205,1e-4);    //from Mathematica
    }

    @Test
    public void testSum() {
        double[] doubleTest = {-1,0,1,2,3};
        long[] longTest = {-1,0,1,2,3};
        byte[] byteTest = {-1,0,1,2,3};
        int[] intTest = {-1,0,1,2,3};
        Assert.assertEquals(MathUtils.sum(doubleTest), 5.0);
        Assert.assertEquals(MathUtils.sum(longTest), 5);
        Assert.assertEquals(MathUtils.sum(byteTest), 5);
        Assert.assertEquals(MathUtils.sum(intTest), 5);
    }

    @Test
    public void testSumRange() {
        long[] longTest = {-1,0,1,2,3};
        Assert.assertEquals(MathUtils.sum(longTest, 1, 4), 3);
    }

    @Test
    public void testMean() {
        double[] test = {0,1,101};
        Assert.assertEquals(MathUtils.mean(test, 0, 0), Double.NaN);
        Assert.assertEquals(MathUtils.mean(test, 0, 1), 0.0);
        Assert.assertEquals(MathUtils.mean(test, 0, 3), 34.0);
    }

    @Test
    public void testStddev() {
        double[] test = {0, -1, 1, -2, 2, -3};
        Assert.assertEquals(MathUtils.stddev(test, 0, 0), Double.NaN);
        Assert.assertEquals(MathUtils.stddev(test, 0, 1), 0.0);
        Assert.assertEquals(MathUtils.stddev(test, 0, 6), 1.707825127659933, 1e-14);
    }

    @Test
    public void testPromote() {
        int[] test = {0, 100, (int)1e10, (int)1e20};
        double[] prom = MathUtils.promote(test);
        for (int i = 0; i < test.length; i++) {
            Assert.assertEquals(test[i], prom[i], 1e-14);
        }
    }

    @Test
    public void testCompareDoubles(){
        Assert.assertEquals(MathUtils.compareDoubles(0.1, 0.2), 1);
        Assert.assertEquals(MathUtils.compareDoubles(0.2, 0.1), -1);
        Assert.assertEquals(MathUtils.compareDoubles(0.1, 0.1), 0);
        Assert.assertEquals(MathUtils.compareDoubles(0.1, 0.1+(1e-7)), 0);
        Assert.assertEquals(MathUtils.compareDoubles(0.1, 0.1+(1e-7), 1e-8), 1);
    }

    @Test
    public void testNormalizeFromReal(){
        final double error = 1e-6;
        final double[] actual = MathUtils.normalizeFromRealSpace(new double[]{1.0, 2.0, 3.0});
        final double[] expected =  {1.0/6.0, 2.0/6.0, 3.0/6.0};
        for (int i = 0; i < actual.length; i++){
            Assert.assertEquals(expected[i], actual[i], error);
        }
    }
}
