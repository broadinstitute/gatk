package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.Arrays;

import static java.lang.Math.log10;

/**
 * Basic unit test for MathUtils
 */
public final class MathUtilsUnitTests {

    private final static Logger logger = LogManager.getLogger(MathUtilsUnitTests.class);

    @BeforeClass
    public void init() {
    }

    @Test
    public void testGoodProbability(){
        Assert.assertTrue(MathUtils.goodLog10Probability(log10(0.1)));
        Assert.assertTrue(MathUtils.goodLog10Probability(log10(0.1), true));
        Assert.assertTrue(MathUtils.goodLog10Probability(log10(0.1), false));
        Assert.assertTrue(MathUtils.goodLog10Probability(Double.NEGATIVE_INFINITY));
        Assert.assertTrue(MathUtils.goodLog10Probability(Double.NEGATIVE_INFINITY, true));
        Assert.assertFalse(MathUtils.goodLog10Probability(Double.NEGATIVE_INFINITY, false));
        Assert.assertTrue(MathUtils.goodProbability(0.1));
        Assert.assertTrue(MathUtils.goodProbability(0.0));
        Assert.assertTrue(MathUtils.goodProbability(1.0));
        Assert.assertFalse(MathUtils.goodProbability(-1.0));
        Assert.assertFalse(MathUtils.goodProbability(2.0));
        Assert.assertFalse(MathUtils.goodProbability(Double.NaN));
        Assert.assertFalse(MathUtils.goodProbability(Double.NEGATIVE_INFINITY));
        Assert.assertFalse(MathUtils.goodProbability(Double.POSITIVE_INFINITY));
    }

    @Test
    public void testLog10OneMinusX(){
        Assert.assertEquals(-0.04575749056, MathUtils.log10OneMinusX(0.1), 1e-6); //result from worlphram alpha
        Assert.assertEquals(0.0, MathUtils.log10OneMinusX(0.0), 1e-6);
        Assert.assertEquals(Double.NEGATIVE_INFINITY, MathUtils.log10OneMinusX(1.0), 1e-6);
    }

    @Test
    public void testLog10Gamma() {
        //The expected values were checked against Wolphram Alpha

        Assert.assertEquals(MathUtils.log10Gamma(4.0), 0.7781513, 1e-6);
        Assert.assertEquals(MathUtils.log10Gamma(10), 5.559763, 1e-6);
        Assert.assertEquals(MathUtils.log10Gamma(10654), 38280.532152137, 1e-6);
    }

    @Test
    public void testBinomialProbability() {
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
    public void testLog10BinomialCoefficient() {
        // note that we can test the binomial coefficient calculation indirectly via Newton's identity
        // (1+z)^m = sum (m choose k)z^k
        double[] z_vals = new double[]{0.999, 0.9, 0.8, 0.5, 0.2, 0.01, 0.0001};
        int[] exponent = new int[]{5, 15, 25, 50, 100};
        for (double z : z_vals) {
            double logz = log10(z);
            for (int exp : exponent) {
                double expected_log = exp * log10(1 + z);
                double[] newtonArray_log = new double[1 + exp];
                for (int k = 0; k <= exp; k++) {
                    newtonArray_log[k] = MathUtils.log10BinomialCoefficient(exp, k) + k * logz;
                }
                Assert.assertEquals(MathUtils.log10sumLog10(newtonArray_log), expected_log, 1e-6);
            }
        }

        Assert.assertEquals(MathUtils.log10BinomialCoefficient(4, 2), 0.7781513, 1e-6);
        Assert.assertEquals(MathUtils.log10BinomialCoefficient(10, 3), 2.079181, 1e-6);
        Assert.assertEquals(MathUtils.log10BinomialCoefficient(103928, 119), 400.2156, 1e-4);
    }


    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLog10BinomialCoefficientErrorN() {
        Assert.assertEquals(MathUtils.log10BinomialCoefficient(-4, 2), 0.7781513, 1e-6);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLog10BinomialCoefficientErrorK() {
        Assert.assertEquals(MathUtils.log10BinomialCoefficient(4, -2), 0.7781513, 1e-6);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testLog10BinomialCoefficientErrorKmoreThanN() {
        Assert.assertEquals(MathUtils.log10BinomialCoefficient(2, 4), 0.7781513, 1e-6);
    }

    @Test
    public void testApproximateLog10SumLog10() {

        final double requiredPrecision = 1E-4;

        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{0.0}), 0.0, requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{-5.15}), -5.15, requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{130.0}), 130.0, requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{-0.145}), -0.145, requiredPrecision);

        Assert.assertEquals(MathUtils.approximateLog10SumLog10(0.0, 0.0), log10(Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-1.0, 0.0), log10(Math.pow(10.0, -1.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(0.0, -1.0), log10(Math.pow(10.0, 0.0) + Math.pow(10.0, -1.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-2.2, -3.5), log10(Math.pow(10.0, -2.2) + Math.pow(10.0, -3.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-1.0, -7.1), log10(Math.pow(10.0, -1.0) + Math.pow(10.0, -7.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(5.0, 6.2), log10(Math.pow(10.0, 5.0) + Math.pow(10.0, 6.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(38.1, 16.2), log10(Math.pow(10.0, 38.1) + Math.pow(10.0, 16.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-38.1, 6.2), log10(Math.pow(10.0, -38.1) + Math.pow(10.0, 6.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-19.1, -37.1), log10(Math.pow(10.0, -19.1) + Math.pow(10.0, -37.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-29.1, -27.6), log10(Math.pow(10.0, -29.1) + Math.pow(10.0, -27.6)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-0.12345, -0.23456), log10(Math.pow(10.0, -0.12345) + Math.pow(10.0, -0.23456)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-15.7654, -17.0101), log10(Math.pow(10.0, -15.7654) + Math.pow(10.0, -17.0101)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-0.12345, Double.NEGATIVE_INFINITY), -0.12345, requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-15.7654, Double.NEGATIVE_INFINITY), -15.7654, requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(Double.NEGATIVE_INFINITY, -15.7654), -15.7654, requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY), Double.NEGATIVE_INFINITY, requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY), Double.POSITIVE_INFINITY, requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY), Double.POSITIVE_INFINITY, requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(Double.NEGATIVE_INFINITY,Double.POSITIVE_INFINITY), Double.POSITIVE_INFINITY, requiredPrecision);

        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{0.0, 0.0}), log10(Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{-1.0, 0.0}), log10(Math.pow(10.0, -1.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{0.0, -1.0}), log10(Math.pow(10.0, 0.0) + Math.pow(10.0, -1.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{-2.2, -3.5}), log10(Math.pow(10.0, -2.2) + Math.pow(10.0, -3.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{-1.0, -7.1}), log10(Math.pow(10.0, -1.0) + Math.pow(10.0, -7.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{5.0, 6.2}), log10(Math.pow(10.0, 5.0) + Math.pow(10.0, 6.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{38.1, 16.2}), log10(Math.pow(10.0, 38.1) + Math.pow(10.0, 16.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{-38.1, 6.2}), log10(Math.pow(10.0, -38.1) + Math.pow(10.0, 6.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{-19.1, -37.1}), log10(Math.pow(10.0, -19.1) + Math.pow(10.0, -37.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{-29.1, -27.6}), log10(Math.pow(10.0, -29.1) + Math.pow(10.0, -27.6)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{-0.12345, -0.23456}), log10(Math.pow(10.0, -0.12345) + Math.pow(10.0, -0.23456)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{-15.7654, -17.0101}), log10(Math.pow(10.0, -15.7654) + Math.pow(10.0, -17.0101)), requiredPrecision);

        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{0.0, 0.0, 0.0}), log10(Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{-1.0, 0.0, 0.0}), log10(Math.pow(10.0, -1.0) + Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{0.0, -1.0, -2.5}), log10(Math.pow(10.0, 0.0) + Math.pow(10.0, -1.0) + Math.pow(10.0, -2.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{-2.2, -3.5, -1.1}), log10(Math.pow(10.0, -2.2) + Math.pow(10.0, -3.5) + Math.pow(10.0, -1.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{-1.0, -7.1, 0.5}), log10(Math.pow(10.0, -1.0) + Math.pow(10.0, -7.1) + Math.pow(10.0, 0.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{5.0, 6.2, 1.3}), log10(Math.pow(10.0, 5.0) + Math.pow(10.0, 6.2) + Math.pow(10.0, 1.3)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{38.1, 16.2, 18.1}), log10(Math.pow(10.0, 38.1) + Math.pow(10.0, 16.2) + Math.pow(10.0, 18.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{-38.1, 6.2, 26.6}), log10(Math.pow(10.0, -38.1) + Math.pow(10.0, 6.2) + Math.pow(10.0, 26.6)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{-19.1, -37.1, -45.1}), log10(Math.pow(10.0, -19.1) + Math.pow(10.0, -37.1) + Math.pow(10.0, -45.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{-29.1, -27.6, -26.2}), log10(Math.pow(10.0, -29.1) + Math.pow(10.0, -27.6) + Math.pow(10.0, -26.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{-0.12345, -0.23456, -0.34567}), log10(Math.pow(10.0, -0.12345) + Math.pow(10.0, -0.23456) + Math.pow(10.0, -0.34567)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{-15.7654, -17.0101, -17.9341}), log10(Math.pow(10.0, -15.7654) + Math.pow(10.0, -17.0101) + Math.pow(10.0, -17.9341)), requiredPrecision);

        Assert.assertEquals(MathUtils.approximateLog10SumLog10(0.0, 0.0, 0.0), log10(Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-1.0, 0.0, 0.0), log10(Math.pow(10.0, -1.0) + Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(0.0, -1.0, -2.5), log10(Math.pow(10.0, 0.0) + Math.pow(10.0, -1.0) + Math.pow(10.0, -2.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-2.2, -3.5, -1.1), log10(Math.pow(10.0, -2.2) + Math.pow(10.0, -3.5) + Math.pow(10.0, -1.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-1.0, -7.1, 0.5), log10(Math.pow(10.0, -1.0) + Math.pow(10.0, -7.1) + Math.pow(10.0, 0.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(5.0, 6.2, 1.3), log10(Math.pow(10.0, 5.0) + Math.pow(10.0, 6.2) + Math.pow(10.0, 1.3)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(38.1, 16.2, 18.1), log10(Math.pow(10.0, 38.1) + Math.pow(10.0, 16.2) + Math.pow(10.0, 18.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-38.1, 6.2, 26.6), log10(Math.pow(10.0, -38.1) + Math.pow(10.0, 6.2) + Math.pow(10.0, 26.6)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-19.1, -37.1, -45.1), log10(Math.pow(10.0, -19.1) + Math.pow(10.0, -37.1) + Math.pow(10.0, -45.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-29.1, -27.6, -26.2), log10(Math.pow(10.0, -29.1) + Math.pow(10.0, -27.6) + Math.pow(10.0, -26.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-0.12345, -0.23456, -0.34567), log10(Math.pow(10.0, -0.12345) + Math.pow(10.0, -0.23456) + Math.pow(10.0, -0.34567)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-15.7654, -17.0101, -17.9341), log10(Math.pow(10.0, -15.7654) + Math.pow(10.0, -17.0101) + Math.pow(10.0, -17.9341)), requiredPrecision);

        // magnitude of the sum doesn't matter, so we can combinatorially test this via partitions of unity
        double[] mult_partitionFactor = new double[]{0.999, 0.98, 0.95, 0.90, 0.8, 0.5, 0.3, 0.1, 0.05, 0.001};
        int[] n_partitions = new int[]{2, 4, 8, 16, 32, 64, 128, 256, 512, 1028};
        for (double alpha : mult_partitionFactor) {
            double log_alpha = log10(alpha);
            double log_oneMinusAlpha = log10(1 - alpha);
            for (int npart : n_partitions) {
                double[] multiplicative = new double[npart];
                double[] equal = new double[npart];
                double remaining_log = 0.0;  // realspace = 1
                for (int i = 0; i < npart - 1; i++) {
                    equal[i] = -log10(npart);
                    double piece = remaining_log + log_alpha; // take a*remaining, leaving remaining-a*remaining = (1-a)*remaining
                    multiplicative[i] = piece;
                    remaining_log = remaining_log + log_oneMinusAlpha;
                }
                equal[npart - 1] = -log10(npart);
                multiplicative[npart - 1] = remaining_log;
                Assert.assertEquals(MathUtils.approximateLog10SumLog10(equal), 0.0, requiredPrecision, String.format("Did not sum to one: k=%d equal partitions.", npart));
                Assert.assertEquals(MathUtils.approximateLog10SumLog10(multiplicative), 0.0, requiredPrecision, String.format("Did not sum to one: k=%d multiplicative partitions with alpha=%f", npart, alpha));
            }
        }
    }

    @Test
    public void testLog10sumLog10() {
        final double requiredPrecision = 1E-14;

        final double log3 = 0.477121254719662;
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{0.0, 0.0, 0.0}), log3, requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{0.0, 0.0, 0.0}, 0), log3, requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{0.0, 0.0, 0.0}, 0, 3), log3, requiredPrecision);

        final double log2 = 0.301029995663981;
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{0.0, 0.0, 0.0}, 0, 2), log2, requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{0.0, 0.0, 0.0}, 0, 1), 0.0, requiredPrecision);

        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{0.0}), 0.0, requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{-5.15}), -5.15, requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{130.0}), 130.0, requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{-0.145}), -0.145, requiredPrecision);

        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{0.0, 0.0}), log10(Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{-1.0, 0.0}), log10(Math.pow(10.0, -1.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{0.0, -1.0}), log10(Math.pow(10.0, 0.0) + Math.pow(10.0, -1.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{-2.2, -3.5}), log10(Math.pow(10.0, -2.2) + Math.pow(10.0, -3.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{-1.0, -7.1}), log10(Math.pow(10.0, -1.0) + Math.pow(10.0, -7.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{5.0, 6.2}), log10(Math.pow(10.0, 5.0) + Math.pow(10.0, 6.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{38.1, 16.2}), log10(Math.pow(10.0, 38.1) + Math.pow(10.0, 16.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{-38.1, 6.2}), log10(Math.pow(10.0, -38.1) + Math.pow(10.0, 6.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{-19.1, -37.1}), log10(Math.pow(10.0, -19.1) + Math.pow(10.0, -37.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{-29.1, -27.6}), log10(Math.pow(10.0, -29.1) + Math.pow(10.0, -27.6)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{-0.12345, -0.23456}), log10(Math.pow(10.0, -0.12345) + Math.pow(10.0, -0.23456)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{-15.7654, -17.0101}), log10(Math.pow(10.0, -15.7654) + Math.pow(10.0, -17.0101)), requiredPrecision);

        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{0.0, 0.0, 0.0}), log10(Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{-1.0, 0.0, 0.0}), log10(Math.pow(10.0, -1.0) + Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{0.0, -1.0, -2.5}), log10(Math.pow(10.0, 0.0) + Math.pow(10.0, -1.0) + Math.pow(10.0, -2.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{-2.2, -3.5, -1.1}), log10(Math.pow(10.0, -2.2) + Math.pow(10.0, -3.5) + Math.pow(10.0, -1.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{-1.0, -7.1, 0.5}), log10(Math.pow(10.0, -1.0) + Math.pow(10.0, -7.1) + Math.pow(10.0, 0.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{5.0, 6.2, 1.3}), log10(Math.pow(10.0, 5.0) + Math.pow(10.0, 6.2) + Math.pow(10.0, 1.3)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{38.1, 16.2, 18.1}), log10(Math.pow(10.0, 38.1) + Math.pow(10.0, 16.2) + Math.pow(10.0, 18.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{-38.1, 6.2, 26.6}), log10(Math.pow(10.0, -38.1) + Math.pow(10.0, 6.2) + Math.pow(10.0, 26.6)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{-19.1, -37.1, -45.1}), log10(Math.pow(10.0, -19.1) + Math.pow(10.0, -37.1) + Math.pow(10.0, -45.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{-29.1, -27.6, -26.2}), log10(Math.pow(10.0, -29.1) + Math.pow(10.0, -27.6) + Math.pow(10.0, -26.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{-0.12345, -0.23456, -0.34567}), log10(Math.pow(10.0, -0.12345) + Math.pow(10.0, -0.23456) + Math.pow(10.0, -0.34567)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{-15.7654, -17.0101, -17.9341}), log10(Math.pow(10.0, -15.7654) + Math.pow(10.0, -17.0101) + Math.pow(10.0, -17.9341)), requiredPrecision);

        // magnitude of the sum doesn't matter, so we can combinatorially test this via partitions of unity
        double[] mult_partitionFactor = new double[]{0.999, 0.98, 0.95, 0.90, 0.8, 0.5, 0.3, 0.1, 0.05, 0.001};
        int[] n_partitions = new int[]{2, 4, 8, 16, 32, 64, 128, 256, 512, 1028};
        for (double alpha : mult_partitionFactor) {
            double log_alpha = log10(alpha);
            double log_oneMinusAlpha = log10(1 - alpha);
            for (int npart : n_partitions) {
                double[] multiplicative = new double[npart];
                double[] equal = new double[npart];
                double remaining_log = 0.0;  // realspace = 1
                for (int i = 0; i < npart - 1; i++) {
                    equal[i] = -log10(npart);
                    double piece = remaining_log + log_alpha; // take a*remaining, leaving remaining-a*remaining = (1-a)*remaining
                    multiplicative[i] = piece;
                    remaining_log = remaining_log + log_oneMinusAlpha;
                }
                equal[npart - 1] = -log10(npart);
                multiplicative[npart - 1] = remaining_log;
                Assert.assertEquals(MathUtils.log10sumLog10(equal), 0.0, requiredPrecision);
                Assert.assertEquals(MathUtils.log10sumLog10(multiplicative), 0.0, requiredPrecision, String.format("Did not sum to one: nPartitions=%d, alpha=%f", npart, alpha));
            }
        }
    }

    @Test
    public void testlog10sumLog10EdgeCases(){
        final double requiredPrecision = 1E-14;

        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{3.0, 2.0, 1.0}, 2, 1), Double.NEGATIVE_INFINITY);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY}), Double.NEGATIVE_INFINITY);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{3.0, 2.0, 1.0}), log10(Math.pow(10.0, 3.0) + Math.pow(10.0, 2.0) + Math.pow(10.0, 1.0)), requiredPrecision);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testlog10sumLog10ErrorNaN(){
        MathUtils.log10sumLog10(new double[]{Double.NaN, 1.0});
    }

    @Test
    public void testNormalize(){
        final double error = 1e-6;
        final double[] normalized = MathUtils.normalizeFromLog10(new double[]{3.0, 2.0, 1.0});
        final double[] normalizedLog = MathUtils.normalizeFromLog10(new double[]{3.0, 2.0, 1.0}, true);
        final double[] normalizedLogInLog = MathUtils.normalizeFromLog10(new double[]{3.0, 2.0, 1.0}, true, true);
        final double[] normalizedExpected    = {0.9009009009009008, 0.09009009009009009, 0.009009009009009009};  //sum of those is 1
        final double[] normalizedLogExpected = {-0.04532297878665748, -1.0453229787866574, -2.0453229787866576}; //sum of 10^ those is 1
        final double[] normalizedLogInLogExpected = {0.0, -1.0, -2.0};
        for (int i = 0; i < normalized.length; i++){
            Assert.assertEquals(normalizedExpected[i], normalized[i], error);
        }
        for (int i = 0; i < normalizedLog.length; i++){
            Assert.assertEquals(normalizedLogExpected[i], normalizedLog[i], error);
        }
        for (int i = 0; i < normalizedLogInLog.length; i++){
            Assert.assertEquals(normalizedLogInLogExpected[i], normalizedLogInLog[i], error);
        }
    }

    @Test
    public void testLog10Factorial() {
        logger.warn("Executing testLog10Factorial");
        Assert.assertEquals(MathUtils.log10Factorial(4), 1.380211, 1e-6);
        Assert.assertEquals(MathUtils.log10Factorial(10), 6.559763, 1e-6);
        Assert.assertEquals(MathUtils.log10Factorial(12), 8.680337, 1e-6);
        Assert.assertEquals(MathUtils.log10Factorial(200), 374.8969, 1e-3);
        Assert.assertEquals(MathUtils.log10Factorial(12342), 45138.26, 1e-1);
        double log10factorial_small = 0;
        double log10factorial_middle = 374.8969;
        double log10factorial_large = 45138.26;
        int small_start = 1;
        int med_start = 200;
        int large_start = 12342;
        for ( int i = 1; i < 1000; i++ ) {
            log10factorial_small += log10(i + small_start);
            log10factorial_middle += log10(i + med_start);
            log10factorial_large += log10(i + large_start);
            Assert.assertEquals(MathUtils.log10Factorial(small_start+i),log10factorial_small,1e-6);
            Assert.assertEquals(MathUtils.log10Factorial(med_start+i),log10factorial_middle,1e-3);
            Assert.assertEquals(MathUtils.log10Factorial(large_start+i),log10factorial_large,1e-1);
        }
    }

    @Test
    public void testCovarianceDivergences() {
        logger.warn("Executing testCovarianceDivergences");
        int size = 3;
        //two symmetric positive-definite matrices
        double[][] cov1 = { {5, 2, 3},
                            {2, 7, 5},
                            {3, 5, 6}};

        double[][] cov2 = { {11, 3, 3},
                            {3, 7, 5},
                            {3, 5, 13}};

        RealMatrix mat1 = new Array2DRowRealMatrix(cov1);
        RealMatrix mat2 = new Array2DRowRealMatrix(cov2);

        Assert.assertEquals(MathUtils.covarianceKLDivergence(mat1, mat2), 3.65393, 1e-4);
        Assert.assertEquals(MathUtils.covarianceGeodesicDistance(mat1, mat2), 1.86205,1e-4);
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
