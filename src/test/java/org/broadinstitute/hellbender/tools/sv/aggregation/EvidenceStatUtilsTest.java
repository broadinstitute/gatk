package org.broadinstitute.hellbender.tools.sv.aggregation;

import com.google.common.collect.Lists;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class EvidenceStatUtilsTest {

    private static final double ERROR_TOL = 1e-6;

    @DataProvider(name = "testProbToQualData")
    public Object[][] testProbToQualData() {
        return new Object[][]{
                {1.0, 99, 1},
                {0., 99, 99},
                {0.1, 99, 10},
                {0.01, 99, 20},
                {0.001, 99, 30},
                {0.001, 15, 15},
                {0.099, 99, 10},
                {0.11, 99, 10}
        };
    }

    @Test(dataProvider= "testProbToQualData")
    public void testProbToQual(final Double pError, final int maxQ, final Integer expected) {
        final Integer test = EvidenceStatUtils.probToQual(pError, (byte) maxQ);
        Assert.assertEquals(test, expected);
    }

    @DataProvider(name = "testCarrierSignalFractionData")
    public Object[][] testCarrierSignalFractionData() {
        return new Object[][]{
                {null, null, null},
                {null, 1., null},
                {1., null, null},
                {0., 0., 0},
                {0., 3., 0},
                {5., 0., 100},
                {4., 4., 50},
                {10., 20., 33},
                {20., 10., 67},
                {1., 49., 2},
                {1., 99., 1},
        };
    }

    @Test(dataProvider= "testCarrierSignalFractionData")
    public void testCarrierSignalFraction(final Double carrierSignal, final Double backgroundSignal, final Integer expected) {
        final Integer test = EvidenceStatUtils.carrierSignalFraction(carrierSignal, backgroundSignal);
        Assert.assertEquals(test, expected);
    }

    @DataProvider(name = "testComputeRepresentativeDepthData")
    public Object[][] testComputeRepresentativeDepthData() {
        return new Object[][]{
                {Lists.newArrayList(10.), 10},
                {Lists.newArrayList(10., 20.), 15},
                {Lists.newArrayList(10., 20.9), 15},
                {Lists.newArrayList(10., 21.), 16},
        };
    }

    @Test(dataProvider= "testComputeRepresentativeDepthData")
    public void testComputeRepresentativeDepth(final Collection<Double> values, int expected) {
        final int test = EvidenceStatUtils.computeRepresentativeDepth(values);
        Assert.assertEquals(test, expected);
    }

    @DataProvider(name = "testCalculateOneSamplePoissonTestData")
    public Object[][] testCalculateOneSamplePoissonTestData() {
        return new Object[][]{
                {new int[]{0}, new double[]{30.}, 0, 1.0},
                {new int[]{0}, new double[]{30.}, 1, 1.0},
                {new int[]{1}, new double[]{30.}, 1, 0.3678794411714422},
                {new int[]{5, 1, 0}, new double[]{30., 35., 40.}, 1, 0.020010047914590948},
        };
    }

    @Test(dataProvider= "testCalculateOneSamplePoissonTestData")
    public void testCalculateOneSamplePoissonTest(final int[] sampleCounts,
                                                  final double[] sampleCoverage,
                                                  final int numCarriers,
                                                  final double expected) {
        final Map<String, Integer> sampleCountsMap = new HashMap<>();
        final Map<String, Double> sampleCoverageMap = new HashMap<>();
        final List<String> carrierSamples = new ArrayList<>(numCarriers);
        final List<String> backgroundSamples = new ArrayList<>(sampleCounts.length - numCarriers);
        for (int i = 0; i < sampleCounts.length; i++) {
            final String sample = "sample" + i;
            sampleCountsMap.put(sample, sampleCounts[i]);
            sampleCoverageMap.put(sample, sampleCoverage[i]);
            // the first numCarriers samples are treated as carriers
            if (i < numCarriers) {
                carrierSamples.add(sample);
            } else {
                backgroundSamples.add(sample);
            }
        }
        final double meanCoverage = MathUtils.sum(sampleCoverage) / sampleCoverage.length;
        final EvidenceStatUtils.PoissonTestResult test = EvidenceStatUtils.calculateOneSamplePoissonTest(sampleCountsMap, carrierSamples,
                backgroundSamples, sampleCoverageMap, meanCoverage);
        Assert.assertTrue(Math.abs(test.getP() - expected) <= ERROR_TOL);
    }

    @DataProvider(name = "testCumulativePoissonProbabilityData")
    public Object[][] testCumulativePoissonProbabilityData() {
        return new Object[][]{
                {1.0, 0},
                {1.0, 1},
                {1.0, 2},
                {1.0, 3},
                {1.0, 10},
                {1.0, 100},
                {2.0, 0},
                {2.0, 1},
                {2.0, 2},
                {2.0, 3},
                {2.0, 10},
                {2.0, 100},
                {2.0, Integer.MAX_VALUE},
                {2.0, -1}
        };
    }

    @Test(dataProvider= "testCumulativePoissonProbabilityData")
    public void testCumulativePoissonProbability(final double lambda, final int val) {
        final double expected = new PoissonDistribution(lambda).cumulativeProbability(val);
        final double test = EvidenceStatUtils.cumulativePoissonProbability(lambda, val);
        Assert.assertTrue(Math.abs(test - expected) <= ERROR_TOL);
    }

    @DataProvider(name = "testGetMedianNormalizedCountData")
    public Object[][] testGetMedianNormalizedCountData() {
        return new Object[][]{
                {new double[]{}, new int[]{}, 0},
                {new double[]{30.}, new int[]{0}, 0},
                {new double[]{30.}, new int[]{1}, 1/30.},
                {new double[]{30., 30.}, new int[]{0, 0}, 0},
                {new double[]{30., 30.}, new int[]{1, 1}, 1/30.},
                {new double[]{30., 45.}, new int[]{1, 1}, 0.5*((1/30.)+(1/45.))},
                {new double[]{30., 30., 30.}, new int[]{1, 0, 2}, 1/30.},
                {new double[]{30., 30., 30.}, new int[]{0, 1, 2}, 1/30.},
                {new double[]{15., 30., 45.}, new int[]{1, 1, 1}, 1/30.},
                {new double[]{16., 30., 45.}, new int[]{1, 2, 2}, 1/16.},
        };
    }

    @Test(dataProvider= "testGetMedianNormalizedCountData")
    public void testGetMedianNormalizedCount(final double[] sampleCoverage, final int[] sampleCounts, final double expected) {
        final Map<String, Double> sampleCoverageMap = new HashMap<>();
        final Map<String, Integer> sampleCountsMap = new HashMap<>();
        final Set<String> samples = new HashSet<>();
        for (int i = 0; i < sampleCoverage.length; i++) {
            final String sample = "sample" + i;
            sampleCoverageMap.put(sample, sampleCoverage[i]);
            if (sampleCounts[i] > 0) {
                sampleCountsMap.put(sample, sampleCounts[i]);
            }
            samples.add(sample);
        }

        final double test = EvidenceStatUtils.getMedianNormalizedCount(samples, sampleCountsMap, sampleCoverageMap);
        Assert.assertTrue(Math.abs(test - expected) <= ERROR_TOL);
    }
}