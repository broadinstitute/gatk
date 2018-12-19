package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class GermlineFilterUnitTest {

    @Test(dataProvider = "log10ProbabilityData")
    public void testGermlineProbability(final double normalLog10Odds, final double log10OddsOfGermlineHetVsSomatic,
                                        final double log10OddsOfGermlineHomAltVsSomatic,
                                        final double populationAF, final double log10SomaticPrior,
                                        final double expectedLog10Posterior, final double tolerance) {
        final double actual = Math.log10(GermlineFilter.germlineProbability(normalLog10Odds, log10OddsOfGermlineHetVsSomatic,
                log10OddsOfGermlineHomAltVsSomatic, populationAF, log10SomaticPrior));
        Assert.assertEquals(actual, expectedLog10Posterior, tolerance);
    }

    @DataProvider(name = "log10ProbabilityData")
    public Object[][] log10ProbabilityData() {
        // normalLog10Odds, log10OddsOfGermlineHetVsSomatic, log10OddsOfGermlineHomAltVsSomatic, populationAF, log10SomaticPrior, expectedLog10Posterior, tolerance
        return new Object[][] {
                // extreme data against normal means not germline, even if the tumor looks like a germline ht i.e. AF = 0.5
                {-100, 0, -10, 0.5, -6, -100, 10},
                // strong evidence in normal, even if population AF is small
                {20, 0, 0, 1e-8, -3, 0, 0.0001},
                //no normal (lod = 0) rare variant
                {0, 0, 0, 1e-10, -6, -4, 2},
                //limit of AF = 1 --> always germline
                {0, -5, -5, 1.0, -3, 0, 0.001},
                //a simple exact value done by hand
                {0, 0, 0, 0.5, -1, -0.01579426, 0.0001},
        };
    }
}