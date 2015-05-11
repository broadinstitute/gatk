package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 * Basic unit test for MathUtils
 */
public final class MathUtilsUnitTests {

    private final static Logger logger = LogManager.getLogger(MathUtilsUnitTests.class);

    @BeforeClass
    public void init() {
    }


    @Test
    public void testLog10Gamma() {
        //The expected values were checked against Wolphram Alpha

        Assert.assertEquals(MathUtils.log10Gamma(4.0), 0.7781513, 1e-6);
        Assert.assertEquals(MathUtils.log10Gamma(10), 5.559763, 1e-6);
        Assert.assertEquals(MathUtils.log10Gamma(10654), 38280.532152137, 1e-6);
    }

    @Test
    public void testLog10BinomialCoefficient() {
        logger.warn("Executing testLog10BinomialCoefficient");
        // note that we can test the binomial coefficient calculation indirectly via Newton's identity
        // (1+z)^m = sum (m choose k)z^k
        double[] z_vals = new double[]{0.999,0.9,0.8,0.5,0.2,0.01,0.0001};
        int[] exponent = new int[]{5,15,25,50,100};
        for ( double z : z_vals ) {
            double logz = Math.log10(z);
            for ( int exp : exponent ) {
                double expected_log = exp*Math.log10(1+z);
                double[] newtonArray_log = new double[1+exp];
                for ( int k = 0 ; k <= exp; k++ ) {
                    newtonArray_log[k] = MathUtils.log10BinomialCoefficient(exp,k)+k*logz;
                }
                Assert.assertEquals(MathUtils.log10sumLog10(newtonArray_log),expected_log,1e-6);
            }
        }

        Assert.assertEquals(MathUtils.log10BinomialCoefficient(4, 2), 0.7781513, 1e-6);
        Assert.assertEquals(MathUtils.log10BinomialCoefficient(10, 3), 2.079181, 1e-6);
        Assert.assertEquals(MathUtils.log10BinomialCoefficient(103928, 119), 400.2156, 1e-4);
    }

    @Test
    public void testLog10sumLog10() {
        final double requiredPrecision = 1E-14;

        final double log3 = 0.477121254719662;
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{0.0, 0.0, 0.0}), log3, requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {0.0, 0.0, 0.0}, 0), log3, requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{0.0, 0.0, 0.0}, 0, 3), log3, requiredPrecision);

        final double log2 = 0.301029995663981;
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {0.0, 0.0, 0.0}, 0, 2), log2, requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {0.0, 0.0, 0.0}, 0, 1), 0.0, requiredPrecision);

        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {0.0}), 0.0, requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-5.15}), -5.15, requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {130.0}), 130.0, requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-0.145}), -0.145, requiredPrecision);

        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {0.0, 0.0}), Math.log10(Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-1.0, 0.0}), Math.log10(Math.pow(10.0, -1.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {0.0, -1.0}), Math.log10(Math.pow(10.0, 0.0) + Math.pow(10.0, -1.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-2.2, -3.5}), Math.log10(Math.pow(10.0, -2.2) + Math.pow(10.0, -3.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-1.0, -7.1}), Math.log10(Math.pow(10.0, -1.0) + Math.pow(10.0, -7.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {5.0, 6.2}), Math.log10(Math.pow(10.0, 5.0) + Math.pow(10.0, 6.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {38.1, 16.2}), Math.log10(Math.pow(10.0, 38.1) + Math.pow(10.0, 16.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-38.1, 6.2}), Math.log10(Math.pow(10.0, -38.1) + Math.pow(10.0, 6.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-19.1, -37.1}), Math.log10(Math.pow(10.0, -19.1) + Math.pow(10.0, -37.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-29.1, -27.6}), Math.log10(Math.pow(10.0, -29.1) + Math.pow(10.0, -27.6)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-0.12345, -0.23456}), Math.log10(Math.pow(10.0, -0.12345) + Math.pow(10.0, -0.23456)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-15.7654, -17.0101}), Math.log10(Math.pow(10.0, -15.7654) + Math.pow(10.0, -17.0101)), requiredPrecision);

        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {0.0, 0.0, 0.0}), Math.log10(Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-1.0, 0.0, 0.0}), Math.log10(Math.pow(10.0, -1.0) + Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {0.0, -1.0, -2.5}), Math.log10(Math.pow(10.0, 0.0) + Math.pow(10.0, -1.0) + Math.pow(10.0, -2.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-2.2, -3.5, -1.1}), Math.log10(Math.pow(10.0, -2.2) + Math.pow(10.0, -3.5) + Math.pow(10.0, -1.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-1.0, -7.1, 0.5}), Math.log10(Math.pow(10.0, -1.0) + Math.pow(10.0, -7.1) + Math.pow(10.0, 0.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {5.0, 6.2, 1.3}), Math.log10(Math.pow(10.0, 5.0) + Math.pow(10.0, 6.2) + Math.pow(10.0, 1.3)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {38.1, 16.2, 18.1}), Math.log10(Math.pow(10.0, 38.1) + Math.pow(10.0, 16.2) + Math.pow(10.0, 18.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-38.1, 6.2, 26.6}), Math.log10(Math.pow(10.0, -38.1) + Math.pow(10.0, 6.2) + Math.pow(10.0, 26.6)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-19.1, -37.1, -45.1}), Math.log10(Math.pow(10.0, -19.1) + Math.pow(10.0, -37.1) + Math.pow(10.0, -45.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-29.1, -27.6, -26.2}), Math.log10(Math.pow(10.0, -29.1) + Math.pow(10.0, -27.6) + Math.pow(10.0, -26.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-0.12345, -0.23456, -0.34567}), Math.log10(Math.pow(10.0, -0.12345) + Math.pow(10.0, -0.23456) + Math.pow(10.0, -0.34567)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-15.7654, -17.0101, -17.9341}), Math.log10(Math.pow(10.0, -15.7654) + Math.pow(10.0, -17.0101) + Math.pow(10.0, -17.9341)), requiredPrecision);

        // magnitude of the sum doesn't matter, so we can combinatorially test this via partitions of unity
        double[] mult_partitionFactor = new double[]{0.999,0.98,0.95,0.90,0.8,0.5,0.3,0.1,0.05,0.001};
        int[] n_partitions = new int[] {2,4,8,16,32,64,128,256,512,1028};
        for ( double alpha : mult_partitionFactor ) {
            double log_alpha = Math.log10(alpha);
            double log_oneMinusAlpha = Math.log10(1-alpha);
            for ( int npart : n_partitions ) {
                double[] multiplicative = new double[npart];
                double[] equal = new double[npart];
                double remaining_log = 0.0;  // realspace = 1
                for ( int i = 0 ; i < npart-1; i++ ) {
                    equal[i] = -Math.log10(npart);
                    double piece = remaining_log + log_alpha; // take a*remaining, leaving remaining-a*remaining = (1-a)*remaining
                    multiplicative[i] = piece;
                    remaining_log = remaining_log + log_oneMinusAlpha;
                }
                equal[npart-1] = -Math.log10(npart);
                multiplicative[npart-1] = remaining_log;
                Assert.assertEquals(MathUtils.log10sumLog10(equal),0.0,requiredPrecision);
                Assert.assertEquals(MathUtils.log10sumLog10(multiplicative),0.0,requiredPrecision,String.format("Did not sum to one: nPartitions=%d, alpha=%f",npart,alpha));
            }
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
            log10factorial_small += Math.log10(i+small_start);
            log10factorial_middle += Math.log10(i+med_start);
            log10factorial_large += Math.log10(i+large_start);
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
}
