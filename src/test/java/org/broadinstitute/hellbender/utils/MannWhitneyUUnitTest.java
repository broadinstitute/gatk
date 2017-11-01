package org.broadinstitute.hellbender.utils;

import com.google.common.primitives.Doubles;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class MannWhitneyUUnitTest extends GATKBaseTest {
    private static double DELTA_PRECISION = 0.00001;

    private static final MannWhitneyU rst = new MannWhitneyU();

    @DataProvider(name="rankSumTestData")
    public Object[][] dataProvider() {
        return new Object[][] {
                new Object[] {"test1", new double[] {20,20,20,20,20}, new double[] {20,20,20,20,21}, 10d},
                new Object[] {"test2", new double[] {20,20,20,20,21}, new double[] {20,20,20,20,21}, 12.5d},
                new Object[] {"test3", new double[] {13,14,15,15,16}, new double[] {16,20,20,21,21}, 0.5d},
                new Object[] {"test4", new double[] {13,14,15,15,16}, new double[] {16,20,20,21,21,21,21,22,23,27}, 0.5d},
                new Object[] {"test5", new double[] {13,14,15,15,16,18,20,22,25,24,25,26,27,28,22,23,19,30,28,22,17},
                        new double[] {16,20,20,21,21,21,21,22,23,27,26,28,29,32,31,22,21,19,16,24,29},
                        180.5d},
                new Object[] {"test6", new double[] {13,14,15,15,16,18,13,14,15,15,16,18,13,14,15,15,16,18,13,14,15,15,16,18},
                        new double[] {21,22,23,27,26,28,29,21,22,23,27,26,28,29,21,22,23,27,26,28,29,21,22,23,27,26,28,29},
                        0d},
                new Object[] {"test7", new double[] {11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20},
                        new double[] {12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20,21},
                        145d},
                new Object[] {"test7", new double[] {11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20},
                        new double[] {12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20,21,21},
                        162d},
                new Object[] {"test8", new double[] {20,20,20,20,20}, new double[] {20,20,20,20,20,20,20,20,20,20}, 25d},
        };
    }

    @Test(dataProvider = "rankSumTestData")
    public void testSimpleU(String name, double[] series1, double[] series2, double U) {
        MannWhitneyU.Result test = rst.test(series1, series2, MannWhitneyU.TestType.TWO_SIDED);
        Assert.assertEquals(test.getU(), U, name);
    }

    @DataProvider(name="oneSidedPTestData")
    public Object[][] oneSidedDataProvider() {
        return new Object[][] {
                new Object[] {"test0", new double[] {0,0}, new double[] {1,1}, 0.083333333},
                new Object[] {"test1", new double[] {20,20,20,20,20}, new double[] {20,20,20,20,21}, .25},
                new Object[] {"test2", new double[] {20,20,20,20,21}, new double[] {20,20,20,20,21}, .5},
                new Object[] {"test3", new double[] {13,14,15,15,16}, new double[] {16,20,20,21,21}, 0.00396825},
                new Object[] {"test4", new double[] {13,14,15,15,16}, new double[] {16,20,20,21,21,21,21,22,23,27}, 0.001469192},
                new Object[] {"test5", new double[] {20,20,20,20,20}, new double[] {20,20,20,20,20,20,20,20,20,20}, .5},
                new Object[] {"test6", new double[] {1,2,3,4,5}, new double[] {6,7,8,9,10}, 0.001984},
                new Object[] {"test7", new double[] {6,7,8,9,10}, new double[] {1,2,3,4,5}, 0.99801587},
                new Object[] {"test8", new double[] {16,20,20,21,21,21,21,22,23,27,16,20,20,21,21,21,21,22,23,27},
                        new double[] {22,23,27,16,20,20,21,21,21,21,22,23,27,60,60}, .08303102},
                new Object[] {"test9", new double[] {16,20,20,21,21,21,21,21,20},
                        new double[] {22,23,27,16,20,20,20,20,21}, 0.388204},
        };
    }

    @Test(dataProvider = "oneSidedPTestData")
    public void testOnesidedP(String name, double[] series1, double[] series2, double P) {
        MannWhitneyU.Result test = rst.test(series1, series2, MannWhitneyU.TestType.FIRST_DOMINATES);
        Assert.assertEquals(test.getP(), P, DELTA_PRECISION, name);
    }

    @DataProvider(name="oneSidedZTestData")
    public Object[][] oneSidedZDataProvider() {
        return new Object[][] {
                new Object[] {"test1", new double[] {20}, new double[] {20,20,20}, 0},
                new Object[] {"test2", new double[] {1}, new double[] {1,2,3}, -0.67448975},
                new Object[] {"test3", new double[] {1,2,3}, new double[] {3}, -0.67448975},
                new Object[] {"test4", new double[] {1,2,3}, new double[] {1}, 0.67448975},
                new Object[] {"test5", new double[] {3,3}, new double[] {1,2,3}, 1.036433},
                new Object[] {"test6", new double[] {20,20,20}, new double[] {20,20,20}, 0},
                new Object[] {"test7", new double[] {20,20,20,20,20}, new double[] {20,20,20,20,20,20,20,20,20,20,20,20,20}, 0},
                new Object[] {"test8", new double[] {1}, new double[] {2}, -0.67448975},
                new Object[] {"test9", new double[] {1}, new double[] {1}, 0},
                new Object[] {"test10", new double[] {60,70,70,60,60,60,60,60}, new double[] {60,60,60,60,60}, .91732119}
        };
    }

    @Test(dataProvider = "oneSidedZTestData")
    public void testOnesidedZ(String name, double[] series1, double[] series2, double Z) {
        MannWhitneyU.Result test = rst.test(series1, series2, MannWhitneyU.TestType.FIRST_DOMINATES);
        Assert.assertEquals(test.getZ(), Z, DELTA_PRECISION, name);
    }

    @Test
    public void testTooManyTies(){
        ArrayList<Integer> listOfNumberOfTies = new ArrayList<>(Arrays.asList(26,3,6,4,13,18,29,36,60,58,87,63,98,125,158,185,193,171,17592,115,100,141,216,298,451,719,1060,1909,3210,5167,7135,10125,11035,3541,732,9));
        Assert.assertEquals(rst.transformTies(64890, listOfNumberOfTies), 8.41378729572e+12);
    }

    @DataProvider(name = "DistributionData")
    public Object[][] makeDistributionData() {
        List<Object[]> tests = new ArrayList<>();

        final int observations = 100;
        final int skew = 3;
        final List<Double> distribution20 = new ArrayList<>(observations);
        final List<Double> distribution30 = new ArrayList<>(observations);
        final List<Double> distribution20_40 = new ArrayList<>(observations);

        makeDistribution(distribution20, 20, skew, observations);
        makeDistribution(distribution30, 30, skew, observations);
        makeDistribution(distribution20_40, 20, skew, observations/2);
        makeDistribution(distribution20_40, 40, skew, observations/2);

        // shuffle the observations
        Utils.resetRandomGenerator();
        Collections.shuffle(distribution20, Utils.getRandomGenerator());
        Collections.shuffle(distribution30, Utils.getRandomGenerator());
        Collections.shuffle(distribution20_40, Utils.getRandomGenerator());

        for ( final int numToReduce : Arrays.asList(0, 10, 50, 100) ) {
            tests.add(new Object[]{distribution20, distribution20, numToReduce, true, observations, "20-20"});
            tests.add(new Object[]{distribution30, distribution30, numToReduce, true, observations, "30-30"});
            tests.add(new Object[]{distribution20_40, distribution20_40, numToReduce, true, observations, "20/40-20/40"});

            tests.add(new Object[]{distribution20, distribution30, numToReduce, false, observations, "20-30"});
            tests.add(new Object[]{distribution30, distribution20, numToReduce, false, observations, "30-20"});

            tests.add(new Object[]{distribution20, distribution20_40, numToReduce, false, observations, "20-20/40"});
            tests.add(new Object[]{distribution30, distribution20_40, numToReduce, true, observations, "30-20/40"});
        }

        return tests.toArray(new Object[][]{});
    }

    private static void makeDistribution(final List<Double> result, final int target, final int skew, final int numObservations) {
        final int rangeStart = target - skew;
        final int rangeEnd = target + skew;

        double current = rangeStart;
        for ( int i = 0; i < numObservations; i++ ) {
            result.add(current++);
            if ( current > rangeEnd )
                current = rangeStart;
        }
    }

    @Test(dataProvider = "DistributionData")
    public void testDistribution(final List<Double> distribution1, final List<Double> distribution2, final int numToReduceIn2, final boolean distributionsShouldBeEqual, final int observations, final String debugString) {
        final MannWhitneyU mannWhitneyU = new MannWhitneyU();

        final List<Double> dist2 = new ArrayList<>(distribution2);
        if ( numToReduceIn2 > 0 ) {
            Double counts = 0.0;
            Double quals = 0.0;

            for ( int i = 0; i < numToReduceIn2; i++ ) {
                counts++;
                quals += dist2.remove(0);
            }

            final Double qual = quals / counts;
            for ( int i = 0; i < numToReduceIn2; i++ )
                dist2.add(qual);
        }

        final Double result = mannWhitneyU.test(Doubles.toArray(distribution1), Doubles.toArray(dist2), MannWhitneyU.TestType.TWO_SIDED).getP();
        Assert.assertFalse(Double.isNaN(result));

        if ( distributionsShouldBeEqual ) {
            if ( numToReduceIn2 >= observations / 2 )
                return;
            Assert.assertTrue(result > 0.1, String.format("%f %d %f", result, numToReduceIn2, dist2.get(0)));
        } else {
            Assert.assertTrue(result < 0.01, String.format("%f %d %f", result, numToReduceIn2, dist2.get(0)));
        }
    }
}
