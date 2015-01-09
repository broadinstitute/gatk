/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.*;

/**
 * Basic unit test for MathUtils
 */
public class MathUtilsUnitTests {

    private final static Logger logger = LogManager.getLogger(MathUtilsUnitTests.class);

    @BeforeClass
    public void init() {
    }

    @Test(dataProvider = "log10OneMinusPow10Data")
    public void testLog10OneMinusPow10(final double x, final double expected) {
        final double actual = MathUtils.log10OneMinusPow10(x);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual,expected,1E-9);
    }

    @Test(dataProvider = "log1mexpData")
    public void testLog1mexp(final double x, final double expected) {
        final double actual = MathUtils.log1mexp(x);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual,expected,1E-9);
    }

    @DataProvider(name = "log10OneMinusPow10Data")
    public Iterator<Object[]> log10OneMinusPow10Data() {

        final double[] inValues = new double[] { Double.NaN, 10, 1, 0, -1, -3, -10, -30, -100, -300, -1000, -3000 };
        return new Iterator<Object[]>() {

            private int i = 0;

            @Override
            public boolean hasNext() {
                return i < inValues.length;

            }

            @Override
            public Object[] next() {
                final double input = inValues[i++];
                final double output = Math.log10( 1 - Math.pow(10,input));
                return new Object[] { input, output };
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @DataProvider(name = "log1mexpData")
    public Iterator<Object[]> log1mexpData() {

        final double[] inValues = new double[] { Double.NaN, 10, 1, 0, -1, -3, -10, -30, -100, -300, -1000, -3000 };
        return new Iterator<Object[]>() {

            private int i = 0;

            @Override
            public boolean hasNext() {
                return i < inValues.length;

            }

            @Override
            public Object[] next() {
                final double input = inValues[i++];
                final double output = Math.log( 1 - Math.exp(input));
                return new Object[] { input, output };
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    /**
     * Tests that we get the right values from the binomial distribution
     */
    @Test
    public void testBinomialProbability() {
        logger.warn("Executing testBinomialProbability");

        Assert.assertEquals(MathUtils.binomialProbability(3, 2, 0.5), 0.375, 0.0001);
        Assert.assertEquals(MathUtils.binomialProbability(100, 10, 0.5), 1.365543e-17, 1e-18);
        Assert.assertEquals(MathUtils.binomialProbability(217, 73, 0.02), 4.521904e-67, 1e-68);
        Assert.assertEquals(MathUtils.binomialProbability(300, 100, 0.02), 9.27097e-91, 1e-92);
        Assert.assertEquals(MathUtils.binomialProbability(300, 150, 0.98), 6.462892e-168, 1e-169);
        Assert.assertEquals(MathUtils.binomialProbability(300, 120, 0.98), 3.090054e-221, 1e-222);
        Assert.assertEquals(MathUtils.binomialProbability(300, 112, 0.98), 2.34763e-236, 1e-237);
    }


    /**
     * Tests that we get the right values from the multinomial distribution
     */
    @Test
    public void testMultinomialProbability() {
        logger.warn("Executing testMultinomialProbability");

        int[] counts0 = {2, 0, 1};
        double[] probs0 = {0.33, 0.33, 0.34};
        Assert.assertEquals(MathUtils.multinomialProbability(counts0, probs0), 0.111078, 1e-6);

        int[] counts1 = {10, 20, 30};
        double[] probs1 = {0.25, 0.25, 0.50};
        Assert.assertEquals(MathUtils.multinomialProbability(counts1, probs1), 0.002870301, 1e-9);

        int[] counts2 = {38, 82, 50, 36};
        double[] probs2 = {0.25, 0.25, 0.25, 0.25};
        Assert.assertEquals(MathUtils.multinomialProbability(counts2, probs2), 1.88221e-09, 1e-10);

        int[] counts3 = {1, 600, 1};
        double[] probs3 = {0.33, 0.33, 0.34};
        Assert.assertEquals(MathUtils.multinomialProbability(counts3, probs3), 5.20988e-285, 1e-286);
    }

    /**
     * Tests that the random index selection is working correctly
     */
    @Test
    public void testRandomIndicesWithReplacement() {
        logger.warn("Executing testRandomIndicesWithReplacement");

        // Check that the size of the list returned is correct
        Assert.assertTrue(MathUtils.sampleIndicesWithReplacement(5, 0).size() == 0);
        Assert.assertTrue(MathUtils.sampleIndicesWithReplacement(5, 1).size() == 1);
        Assert.assertTrue(MathUtils.sampleIndicesWithReplacement(5, 5).size() == 5);
        Assert.assertTrue(MathUtils.sampleIndicesWithReplacement(5, 1000).size() == 1000);

        // Check that the list contains only the k element range that as asked for - no more, no less
        List<Integer> Five = new ArrayList<>();
        Collections.addAll(Five, 0, 1, 2, 3, 4);
        List<Integer> BigFive = MathUtils.sampleIndicesWithReplacement(5, 10000);
        Assert.assertTrue(BigFive.containsAll(Five));
        Assert.assertTrue(Five.containsAll(BigFive));
    }

    /**
     * Tests that we get the right values from the multinomial distribution
     */
    @Test
    public void testSliceListByIndices() {
        logger.warn("Executing testSliceListByIndices");

        // Check that the list contains only the k element range that as asked for - no more, no less but now
        // use the index list to pull elements from another list using sliceListByIndices
        List<Integer> Five = new ArrayList<>();
        Collections.addAll(Five, 0, 1, 2, 3, 4);
        List<Character> FiveAlpha = new ArrayList<>();
        Collections.addAll(FiveAlpha, 'a', 'b', 'c', 'd', 'e');
        List<Integer> BigFive = MathUtils.sampleIndicesWithReplacement(5, 10000);
        List<Character> BigFiveAlpha = MathUtils.sliceListByIndices(BigFive, FiveAlpha);
        Assert.assertTrue(BigFiveAlpha.containsAll(FiveAlpha));
        Assert.assertTrue(FiveAlpha.containsAll(BigFiveAlpha));
    }

    /**
     * Tests that we correctly compute mean and standard deviation from a stream of numbers
     */
    @Test
    public void testRunningAverage() {
        logger.warn("Executing testRunningAverage");

        int[] numbers = {1, 2, 4, 5, 3, 128, 25678, -24};
        MathUtils.RunningAverage r = new MathUtils.RunningAverage();

        for (final double b : numbers)
            r.add(b);

        Assert.assertEquals((long) numbers.length, r.observationCount());
        Assert.assertTrue(r.mean() - 3224.625 < 2e-10);
        Assert.assertTrue(r.stddev() - 9072.6515881128 < 2e-10);
    }

    @Test
    public void testLog10Gamma() {
        logger.warn("Executing testLog10Gamma");

        Assert.assertEquals(MathUtils.log10Gamma(4.0), 0.7781513, 1e-6);
        Assert.assertEquals(MathUtils.log10Gamma(10), 5.559763, 1e-6);
        Assert.assertEquals(MathUtils.log10Gamma(10654), 38280.53, 1e-2);
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
    public void testFactorial() {
        logger.warn("Executing testFactorial");
        Assert.assertEquals((int) MathUtils.factorial(4), 24);
        Assert.assertEquals((int) MathUtils.factorial(10), 3628800);
        Assert.assertEquals((int) MathUtils.factorial(12), 479001600);
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
    public void testApproximateLog10SumLog10() {

        final double requiredPrecision = 1E-4;

        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {0.0}), 0.0, requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {-5.15}), -5.15, requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {130.0}), 130.0, requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {-0.145}), -0.145, requiredPrecision);

        Assert.assertEquals(MathUtils.approximateLog10SumLog10(0.0, 0.0), Math.log10(Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-1.0, 0.0), Math.log10(Math.pow(10.0, -1.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(0.0, -1.0), Math.log10(Math.pow(10.0, 0.0) + Math.pow(10.0, -1.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-2.2, -3.5), Math.log10(Math.pow(10.0, -2.2) + Math.pow(10.0, -3.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-1.0, -7.1), Math.log10(Math.pow(10.0, -1.0) + Math.pow(10.0, -7.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(5.0, 6.2), Math.log10(Math.pow(10.0, 5.0) + Math.pow(10.0, 6.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(38.1, 16.2), Math.log10(Math.pow(10.0, 38.1) + Math.pow(10.0, 16.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-38.1, 6.2), Math.log10(Math.pow(10.0, -38.1) + Math.pow(10.0, 6.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-19.1, -37.1), Math.log10(Math.pow(10.0, -19.1) + Math.pow(10.0, -37.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-29.1, -27.6), Math.log10(Math.pow(10.0, -29.1) + Math.pow(10.0, -27.6)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-0.12345, -0.23456), Math.log10(Math.pow(10.0, -0.12345) + Math.pow(10.0, -0.23456)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-15.7654, -17.0101), Math.log10(Math.pow(10.0, -15.7654) + Math.pow(10.0, -17.0101)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-0.12345, Double.NEGATIVE_INFINITY), -0.12345, requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-15.7654, Double.NEGATIVE_INFINITY), -15.7654, requiredPrecision);

        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {0.0, 0.0}), Math.log10(Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {-1.0, 0.0}), Math.log10(Math.pow(10.0, -1.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {0.0, -1.0}), Math.log10(Math.pow(10.0, 0.0) + Math.pow(10.0, -1.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {-2.2, -3.5}), Math.log10(Math.pow(10.0, -2.2) + Math.pow(10.0, -3.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {-1.0, -7.1}), Math.log10(Math.pow(10.0, -1.0) + Math.pow(10.0, -7.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {5.0, 6.2}), Math.log10(Math.pow(10.0, 5.0) + Math.pow(10.0, 6.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {38.1, 16.2}), Math.log10(Math.pow(10.0, 38.1) + Math.pow(10.0, 16.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {-38.1, 6.2}), Math.log10(Math.pow(10.0, -38.1) + Math.pow(10.0, 6.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {-19.1, -37.1}), Math.log10(Math.pow(10.0, -19.1) + Math.pow(10.0, -37.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {-29.1, -27.6}), Math.log10(Math.pow(10.0, -29.1) + Math.pow(10.0, -27.6)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {-0.12345, -0.23456}), Math.log10(Math.pow(10.0, -0.12345) + Math.pow(10.0, -0.23456)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {-15.7654, -17.0101}), Math.log10(Math.pow(10.0, -15.7654) + Math.pow(10.0, -17.0101)), requiredPrecision);

        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {0.0, 0.0, 0.0}), Math.log10(Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {-1.0, 0.0, 0.0}), Math.log10(Math.pow(10.0, -1.0) + Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {0.0, -1.0, -2.5}), Math.log10(Math.pow(10.0, 0.0) + Math.pow(10.0, -1.0) + Math.pow(10.0, -2.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {-2.2, -3.5, -1.1}), Math.log10(Math.pow(10.0, -2.2) + Math.pow(10.0, -3.5) + Math.pow(10.0, -1.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {-1.0, -7.1, 0.5}), Math.log10(Math.pow(10.0, -1.0) + Math.pow(10.0, -7.1) + Math.pow(10.0, 0.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {5.0, 6.2, 1.3}), Math.log10(Math.pow(10.0, 5.0) + Math.pow(10.0, 6.2) + Math.pow(10.0, 1.3)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {38.1, 16.2, 18.1}), Math.log10(Math.pow(10.0, 38.1) + Math.pow(10.0, 16.2) + Math.pow(10.0, 18.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {-38.1, 6.2, 26.6}), Math.log10(Math.pow(10.0, -38.1) + Math.pow(10.0, 6.2) + Math.pow(10.0, 26.6)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {-19.1, -37.1, -45.1}), Math.log10(Math.pow(10.0, -19.1) + Math.pow(10.0, -37.1) + Math.pow(10.0, -45.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {-29.1, -27.6, -26.2}), Math.log10(Math.pow(10.0, -29.1) + Math.pow(10.0, -27.6) + Math.pow(10.0, -26.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {-0.12345, -0.23456, -0.34567}), Math.log10(Math.pow(10.0, -0.12345) + Math.pow(10.0, -0.23456) + Math.pow(10.0, -0.34567)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[] {-15.7654, -17.0101, -17.9341}), Math.log10(Math.pow(10.0, -15.7654) + Math.pow(10.0, -17.0101) + Math.pow(10.0, -17.9341)), requiredPrecision);

        Assert.assertEquals(MathUtils.approximateLog10SumLog10(0.0, 0.0, 0.0), Math.log10(Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-1.0, 0.0, 0.0), Math.log10(Math.pow(10.0, -1.0) + Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(0.0, -1.0, -2.5), Math.log10(Math.pow(10.0, 0.0) + Math.pow(10.0, -1.0) + Math.pow(10.0, -2.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-2.2, -3.5, -1.1), Math.log10(Math.pow(10.0, -2.2) + Math.pow(10.0, -3.5) + Math.pow(10.0, -1.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-1.0, -7.1, 0.5), Math.log10(Math.pow(10.0, -1.0) + Math.pow(10.0, -7.1) + Math.pow(10.0, 0.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(5.0, 6.2, 1.3), Math.log10(Math.pow(10.0, 5.0) + Math.pow(10.0, 6.2) + Math.pow(10.0, 1.3)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(38.1, 16.2, 18.1), Math.log10(Math.pow(10.0, 38.1) + Math.pow(10.0, 16.2) + Math.pow(10.0, 18.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-38.1, 6.2, 26.6), Math.log10(Math.pow(10.0, -38.1) + Math.pow(10.0, 6.2) + Math.pow(10.0, 26.6)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-19.1, -37.1, -45.1), Math.log10(Math.pow(10.0, -19.1) + Math.pow(10.0, -37.1) + Math.pow(10.0, -45.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-29.1, -27.6, -26.2), Math.log10(Math.pow(10.0, -29.1) + Math.pow(10.0, -27.6) + Math.pow(10.0, -26.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-0.12345, -0.23456, -0.34567), Math.log10(Math.pow(10.0, -0.12345) + Math.pow(10.0, -0.23456) + Math.pow(10.0, -0.34567)), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(-15.7654, -17.0101, -17.9341), Math.log10(Math.pow(10.0, -15.7654) + Math.pow(10.0, -17.0101) + Math.pow(10.0, -17.9341)), requiredPrecision);

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
                Assert.assertEquals(MathUtils.approximateLog10SumLog10(equal),0.0,requiredPrecision,String.format("Did not sum to one: k=%d equal partitions.",npart));
                Assert.assertEquals(MathUtils.approximateLog10SumLog10(multiplicative),0.0,requiredPrecision, String.format("Did not sum to one: k=%d multiplicative partitions with alpha=%f",npart,alpha));
            }
        }
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
    public void testLogDotProduct() {
        Assert.assertEquals(MathUtils.logDotProduct(new double[]{-5.0,-3.0,2.0}, new double[]{6.0,7.0,8.0}),10.0,1e-3);
        Assert.assertEquals(MathUtils.logDotProduct(new double[]{-5.0}, new double[]{6.0}),1.0,1e-3);
    }

    @Test
    public void testNormalDistribution() {
        final double requiredPrecision = 1E-10;

        for( final double mu : new double[]{-5.0, -3.2, -1.5, 0.0, 1.2, 3.0, 5.8977} ) {
            for( final double sigma : new double[]{1.2, 3.0, 5.8977} ) {
                for( final double x : new double[]{-5.0, -3.2, -1.5, 0.0, 1.2, 3.0, 5.8977} ) {
                    NormalDistribution n = new NormalDistribution(mu, sigma);
                    Assert.assertEquals(n.density(x), MathUtils.normalDistribution(mu, sigma, x), requiredPrecision, "mu:" + mu + " sigma:" + sigma + " x:" + x);
                    Assert.assertEquals(Math.log10(n.density(x)), MathUtils.normalDistributionLog10(mu, sigma, x), requiredPrecision, "mu:" + mu + " sigma:" + sigma + " x:" + x);
                }
            }
        }
    }

    @DataProvider(name = "ArrayMinData")
    public Object[][] makeArrayMinData() {
        List<Object[]> tests = new ArrayList<>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        tests.add(new Object[]{Arrays.asList(10), 10});
        tests.add(new Object[]{Arrays.asList(-10), -10});

        for ( final List<Integer> values : Utils.makePermutations(Arrays.asList(1,2,3), 3, false) ) {
            tests.add(new Object[]{values, 1});
        }

        for ( final List<Integer> values : Utils.makePermutations(Arrays.asList(1,2,-3), 3, false) ) {
            tests.add(new Object[]{values, -3});
        }


        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ArrayMinData")
    public void testArrayMinList(final List<Integer> values, final int expected) {
        final int actual = MathUtils.arrayMin(values);
        Assert.assertEquals(actual, expected, "Failed with " + values);
    }

    @Test(dataProvider = "ArrayMinData")
    public void testArrayMinIntArray(final List<Integer> values, final int expected) {
        final int[] asArray = ArrayUtils.toPrimitive(values.toArray(new Integer[values.size()]));
        final int actual = MathUtils.arrayMin(asArray);
        Assert.assertEquals(actual, expected, "Failed with " + values);
    }

    @Test(dataProvider = "ArrayMinData")
    public void testArrayMinByteArray(final List<Integer> values, final int expected) {
        final byte[] asArray = new byte[values.size()];
        for ( int i = 0; i < values.size(); i++ ) asArray[i] = (byte)(values.get(i) & 0xFF);
        final byte actual = MathUtils.arrayMin(asArray);
        Assert.assertEquals(actual, (byte)(expected & 0xFF), "Failed with " + values);
    }

    @Test(dataProvider = "ArrayMinData")
    public void testArrayMinDoubleArray(final List<Integer> values, final int expected) {
        final double[] asArray = new double[values.size()];
        for ( int i = 0; i < values.size(); i++ ) asArray[i] = (double)(values.get(i));
        final double actual = MathUtils.arrayMin(asArray);
        Assert.assertEquals(actual, (double)expected, "Failed with " + values);
    }

    @DataProvider(name = "MedianData")
    public Object[][] makeMedianData() {
        final List<Object[]> tests = new ArrayList<>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        tests.add(new Object[]{Arrays.asList(10), 10});
        tests.add(new Object[]{Arrays.asList(1, 10), 10});

        for ( final List<Integer> values : Utils.makePermutations(Arrays.asList(1,2,-3), 3, false) ) {
            tests.add(new Object[]{values, 1});
        }

        for ( final List<Double> values : Utils.makePermutations(Arrays.asList(1.1,2.1,-3.1), 3, false) ) {
            tests.add(new Object[]{values, 1.1});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MedianData")
    @SuppressWarnings({"rawtypes", "unchecked"})
    public void testMedian(final List<Comparable> values, final Comparable expected) {
        final Comparable actual = MathUtils.median(values);
        Assert.assertEquals(actual, expected, "Failed with " + values);
    }



    // man. All this to test dirichlet.

    private double[] unwrap(List<Double> stuff) {
        double[] unwrapped = new double[stuff.size()];
        int idx = 0;
        for ( Double d : stuff ) {
            unwrapped[idx++] = d == null ? 0.0 : d;
        }

        return unwrapped;
    }

    /**
     * The PartitionGenerator generates all of the partitions of a number n, e.g.
     * 5 + 0
     * 4 + 1
     * 3 + 2
     * 3 + 1 + 1
     * 2 + 2 + 1
     * 2 + 1 + 1 + 1
     * 1 + 1 + 1 + 1 + 1
     *
     * This is used to help enumerate the state space over which the Dirichlet-Multinomial is defined,
     * to ensure that the distribution function is properly implemented
     */
    class PartitionGenerator implements Iterator<List<Integer>> {
        // generate the partitions of an integer, each partition sorted numerically
        int n;
        List<Integer> a;

        int y;
        int k;
        int state;

        int x;
        int l;

        public PartitionGenerator(int n) {
            this.n = n;
            this.y = n - 1;
            this.k = 1;
            this.a = new ArrayList<>();
            for ( int i = 0; i < n; i++ ) {
                this.a.add(i);
            }
            this.state = 0;
        }

        public void remove()  { /* do nothing */ }

        public boolean hasNext() { return ! ( this.k == 0 && state == 0 ); }

        private String dataStr()  {
            return String.format("a = [%s]  k = %d  y = %d  state = %d  x = %d  l = %d",
                    Utils.join(",",a), k, y, state, x, l);
        }

        public List<Integer> next() {
            if ( this.state == 0 ) {
                this.x = a.get(k-1)+1;
                k -= 1;
                this.state = 1;
            }

            if ( this.state == 1 ) {
                while ( 2 * x <= y ) {
                    this.a.set(k,x);
                    this.y -= x;
                    this.k++;
                }
                this.l = 1+this.k;
                this.state = 2;
            }

            if ( this.state == 2 ) {
                if ( x <= y ) {
                    this.a.set(k,x);
                    this.a.set(l,y);
                    x += 1;
                    y -= 1;
                    return this.a.subList(0, this.k + 2);
                } else {
                    this.state =3;
                }
            }

            if ( this.state == 3 ) {
                this.a.set(k,x+y);
                this.y = x + y - 1;
                this.state = 0;
                return a.subList(0, k + 1);
            }

            throw new IllegalStateException("Cannot get here");
        }

        public String toString() {
            final StringBuilder buf = new StringBuilder();
            buf.append("{ ");
            while ( hasNext() ) {
                buf.append("[");
                buf.append(Utils.join(",",next()));
                buf.append("],");
            }
            buf.deleteCharAt(buf.lastIndexOf(","));
            buf.append(" }");
            return buf.toString();
        }

    }

    /**
     * NextCounts is the enumerator over the state space of the multinomial dirichlet.
     *
     * It filters the partition of the total sum to only those with a number of terms
     * equal to the number of categories.
     *
     * It then generates all permutations of that partition.
     *
     * In so doing it enumerates over the full state space.
     */
    class NextCounts implements Iterator<int[]> {

        private PartitionGenerator partitioner;
        private int numCategories;
        private int[] next;

        public NextCounts(int numCategories, int totalCounts) {
            partitioner = new PartitionGenerator(totalCounts);
            this.numCategories = numCategories;
            next = nextFromPartitioner();
        }

        public void remove() { /* do nothing */ }

        public boolean hasNext() { return next != null; }

        public int[] next() {
            int[] toReturn = clone(next);
            next = nextPermutation();
            if ( next == null ) {
                next = nextFromPartitioner();
            }

            return toReturn;
        }

        private int[] clone(int[] arr) {
            return Arrays.copyOf(arr, arr.length);
        }

        private int[] nextFromPartitioner() {
            if ( partitioner.hasNext() ) {
                List<Integer> nxt = partitioner.next();
                while ( partitioner.hasNext() && nxt.size() > numCategories ) {
                    nxt = partitioner.next();
                }

                if ( nxt.size() > numCategories ) {
                    return null;
                } else {
                    int[] buf = new int[numCategories];
                    for ( int idx = 0; idx < nxt.size(); idx++ ) {
                        buf[idx] = nxt.get(idx);
                    }
                    Arrays.sort(buf);
                    return buf;
                }
            }

            return null;
        }

        public int[] nextPermutation() {
            return MathUtilsUnitTests.nextPermutation(next);
        }

    }

    public static int[] nextPermutation(int[] next) {
        // the counts can swap among each other. The int[] is originally in ascending order
        // this generates the next array in lexicographic order descending

        // locate the last occurrence where next[k] < next[k+1]
        int gt = -1;
        for ( int idx = 0; idx < next.length-1; idx++) {
            if ( next[idx] < next[idx+1] ) {
                gt = idx;
            }
        }

        if ( gt == -1 ) {
            return null;
        }

        int largestLessThan = gt+1;
        for ( int idx = 1 + largestLessThan; idx < next.length; idx++) {
            if ( next[gt] < next[idx] ) {
                largestLessThan = idx;
            }
        }

        int val = next[gt];
        next[gt] = next[largestLessThan];
        next[largestLessThan] = val;

        // reverse the tail of the array
        int[] newTail = new int[next.length-gt-1];
        int ctr = 0;
        for ( int idx = next.length-1; idx > gt; idx-- ) {
            newTail[ctr++] = next[idx];
        }

        for ( int idx = 0; idx < newTail.length; idx++) {
            next[gt+idx+1] = newTail[idx];
        }

        return next;
    }


    // before testing the dirichlet multinomial, we need to test the
    // classes used to test the dirichlet multinomial

    @Test
    public void testPartitioner() {
        int[] numsToTest = new int[]{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
        int[] expectedSizes = new int[]{1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231, 297, 385, 490, 627};
        for ( int testNum = 0; testNum < numsToTest.length; testNum++ ) {
            PartitionGenerator gen = new PartitionGenerator(numsToTest[testNum]);
            int size = 0;
            while ( gen.hasNext() ) {
                logger.debug(gen.dataStr());
                size += 1;
                gen.next();
            }
            Assert.assertEquals(size,expectedSizes[testNum],
                    String.format("Expected %d partitions, observed %s",expectedSizes[testNum],new PartitionGenerator(numsToTest[testNum]).toString()));
        }
    }

    @Test
    public void testNextPermutation() {
        int[] arr = new int[]{1,2,3,4};
        int[][] gens = new int[][] {
                new int[]{1,2,3,4},
                new int[]{1,2,4,3},
                new int[]{1,3,2,4},
                new int[]{1,3,4,2},
                new int[]{1,4,2,3},
                new int[]{1,4,3,2},
                new int[]{2,1,3,4},
                new int[]{2,1,4,3},
                new int[]{2,3,1,4},
                new int[]{2,3,4,1},
                new int[]{2,4,1,3},
                new int[]{2,4,3,1},
                new int[]{3,1,2,4},
                new int[]{3,1,4,2},
                new int[]{3,2,1,4},
                new int[]{3,2,4,1},
                new int[]{3,4,1,2},
                new int[]{3,4,2,1},
                new int[]{4,1,2,3},
                new int[]{4,1,3,2},
                new int[]{4,2,1,3},
                new int[]{4,2,3,1},
                new int[]{4,3,1,2},
                new int[]{4,3,2,1} };
        for ( int gen = 0; gen < gens.length; gen ++ ) {
            for ( int idx = 0; idx < 3; idx++ ) {
                Assert.assertEquals(arr[idx],gens[gen][idx],
                        String.format("Error at generation %d, expected %s, observed %s",gen,Arrays.toString(gens[gen]),Arrays.toString(arr)));
            }
            arr = nextPermutation(arr);
        }
    }

    private double[] addEpsilon(double[] counts) {
        double[] d = new double[counts.length];
        for ( int i = 0; i < counts.length; i ++ ) {
            d[i] = counts[i] + 1e-3;
        }
        return d;
    }

    @Test
    public void testDirichletMultinomial() {
        List<double[]> testAlleles = Arrays.asList(
                new double[]{80,240},
                new double[]{1,10000},
                new double[]{0,500},
                new double[]{5140,20480},
                new double[]{5000,800,200},
                new double[]{6,3,1000},
                new double[]{100,400,300,800},
                new double[]{8000,100,20,80,2},
                new double[]{90,20000,400,20,4,1280,720,1}
        );

        Assert.assertTrue(! Double.isInfinite(MathUtils.log10Gamma(1e-3)) && ! Double.isNaN(MathUtils.log10Gamma(1e-3)));

        int[] numAlleleSampled = new int[]{2,5,10,20,25};
        for ( double[] alleles : testAlleles ) {
            for ( int count : numAlleleSampled ) {
                // test that everything sums to one. Generate all multinomial draws
                List<Double> likelihoods = new ArrayList<>(100000);
                NextCounts generator = new NextCounts(alleles.length,count);
                double maxLog = Double.MIN_VALUE;
                //List<String> countLog = new ArrayList<String>(200);
                while ( generator.hasNext() ) {
                    int[] thisCount = generator.next();
                    //countLog.add(Arrays.toString(thisCount));
                    Double likelihood = MathUtils.dirichletMultinomial(addEpsilon(alleles),thisCount);
                    Assert.assertTrue(! Double.isNaN(likelihood) && ! Double.isInfinite(likelihood),
                            String.format("Likelihood for counts %s and nAlleles %d was %s",
                                    Arrays.toString(thisCount),alleles.length,Double.toString(likelihood)));
                    if ( likelihood > maxLog )
                        maxLog = likelihood;
                    likelihoods.add(likelihood);
                }
                //System.out.printf("%d likelihoods and max is (probability) %e\n",likelihoods.size(),Math.pow(10,maxLog));
                Assert.assertEquals(MathUtils.sumLog10(unwrap(likelihoods)),1.0,1e-7,
                        String.format("Counts %d and alleles %d have nLikelihoods %d. \n Counts: %s",
                                count,alleles.length,likelihoods.size(), "NODEBUG"/*,countLog*/));
            }
        }
    }
}
