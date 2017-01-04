package org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class ArtifactStatisticsScorerUnitTest extends BaseTest {
    @Test(dataProvider = "BasicTest")
    public void testCalculateSuppressionFactorFromPreAdapterQ(double biasP1, double biasP2, double preAdapterQ, double gt) {
        Assert.assertEquals(ArtifactStatisticsScorer.calculateSuppressionFactorFromPreAdapterQ(preAdapterQ, biasP1, biasP2),gt, 1e-6);
    }

    @DataProvider(name = "BasicTest")
    public Object[][] basicArtifactModesWithRC() {
        return new Object[][]{
                {30.0, 1.5, 20, 0.999999694097773},
                {30.0, 1.5, 21, 0.999998629042793},
                {30.0, 1.5, 22, 0.999993855825398},
                {30.0, 1.5, 23, 0.999972464308885},
                {30.0, 1.5, 24, 0.999876605424014},
                {30.0, 1.5, 25, 0.999447221363076},
                {30.0, 1.5, 26, 0.997527376843365},
                {30.0, 1.5, 27, 0.989013057369407},
                {30.0, 1.5, 28, 0.952574126822433},
                {30.0, 1.5, 29, 0.817574476193644},
                {30.0, 1.5, 30, 0.500000000000000},
                {30.0, 1.5, 31, 0.182425523806356},
                {30.0, 1.5, 32, 0.0474258731775668},
                {30.0, 1.5, 33, 0.0109869426305932},
                {30.0, 1.5, 34, 0.00247262315663477},
                {30.0, 1.5, 35, 0.000552778636923600},
                {30.0, 1.5, 36, 0.000123394575986232},
                {30.0, 1.5, 37, 2.75356911145835e-05},
                {30.0, 1.5, 38, 6.14417460221472e-06},
                {30.0, 1.5, 39, 1.37095720685784e-06},
                {30.0, 1.5, 40, 3.05902226925625e-07},

                // if biasP2 = 0, then result is always 1.0
                // biasQP2=0 forces a sharp cutoff for filtering only OxoG<biasQP1 samples.
                {30.0, 0, 10, 1.0},
                {30.0, 0, 29.9, 1.0},
                {30.0, 0, 30.1, 0.0},
                {30.0, 0, 40, 0.0},
                {30.0, 0, 100, 0.0},
                {50.0, 0, 10, 1.0},
                {50.0, 0, 49.9, 1.0},
                {50.0, 0, 50.1, 0.0},
                {50.0, 0, 40, 1.0},
                {50.0, 0, 100, 0.0},

                // Different values
                {40.0, 1.0, 25, 0.999999694097773},
                {40.0, 1.0, 26, 0.999999168471972},
                {40.0, 1.0, 27, 0.999997739675702},
                {40.0, 1.0, 28, 0.999993855825398},
                {40.0, 1.0, 29, 0.999983298578152},
                {40.0, 1.0, 30, 0.999954602131298},
                {40.0, 1.0, 31, 0.999876605424014},
                {40.0, 1.0, 32, 0.999664649869534},
                {40.0, 1.0, 33, 0.999088948805599},
                {40.0, 1.0, 34, 0.997527376843365},
                {40.0, 1.0, 35, 0.993307149075715}
        };
    }

}
