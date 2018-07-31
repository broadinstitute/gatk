package org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class BetaBinomialDistributionUnitTest extends GATKBaseTest {

    @Test(dataProvider = "basicBetaBinomialTest")
    public void basicTestPdf(int x, int n, int a, int b, double gt) {
        double test = new BetaBinomialDistribution(null, a, b, n).probability(x);
        Assert.assertEquals(test, gt, 1e-3);
    }

    @DataProvider(name = "basicBetaBinomialTest")
    public Object[][] createBasicBetaBinomialTest() {
        return new Object[][] {
                {0, 25, 5, 36, 0.0797},
                {1, 25, 5, 36, 0.1660},
                {2, 25, 5, 36, 0.2025},
                {3, 25, 5, 36, 0.1874},
                {4, 25, 5, 36, 0.1447},
                {5, 25, 5, 36, 0.0976},
                {6, 25, 5, 36, 0.0592},
                {7, 25, 5, 36, 0.0327},
                {8, 25, 5, 36, 0.0167},
                {9, 25, 5, 36, 0.0079},
                {10, 25, 5, 36, 0.0035},
                {11, 25, 5, 36, 0.0014},
                {12, 25, 5, 36, 0.0005},
                {13, 25, 5, 36, 0.0002},
                {14, 25, 5, 36, 0.0001},
                {15, 25, 5, 36, 0.0000},
                {16, 25, 5, 36, 0.0000},
                {17, 25, 5, 36, 0.0000},
                {18, 25, 5, 36, 0.0000},
                {19, 25, 5, 36, 0.0000},
                {20, 25, 5, 36, 0.0000}
        };
    }
}
