package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.broadinstitute.hellbender.tools.exome.AllelicCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Test the initialization on the allele fraction model
 * @author David Benjamin
 */
public final class AlleleFractionInitializerUnitTest {
    @Test
    public void testInitialize() {
        final double averageHetsPerSegment = 20;
        final int numSegments = 100;
        final int averageDepth = 50;

        final double biasMean = 1.1;
        final double biasVariance = 0.01;
        final double outlierProbability = 0.02;

        final double minorFractionTolerance = 0.02;
        final double biasMeanTolerance = 0.02;
        final double biasVarianceTolerance = 0.05;
        //the likelihood is very insensitive to this parameter, so (1) it's posterior is very broad and
        // (2) it's not important to initialize it well (or to sample it well, for that matter)
        // additionally, unlike the other parameters, which converge to the posterior mode more more quickly via initialization
        // than via sampling, this parameter gets updated very efficiently from sampling and slowly in initialization
        // therefore, we really don't care how well we initialize it as long as we do a good job with the other parameters.

        final double outlierProbabilityTolerance = 0.03;
        final AlleleFractionSimulatedData simulatedData = new AlleleFractionSimulatedData(averageHetsPerSegment, numSegments,
                averageDepth, biasMean, biasVariance, outlierProbability);

        final AlleleFractionData data = new AlleleFractionData(simulatedData.getSegmentedModel());
        final AlleleFractionState initializedState = new AlleleFractionInitializer(data).getInitializedState();

        final AlleleFractionSimulatedData.AlleleFractionStateError error = simulatedData.error(initializedState);
        Assert.assertEquals(error.averageMinorFractionError, 0, minorFractionTolerance);
        Assert.assertEquals(error.biasMeanError, 0, biasMeanTolerance);
        Assert.assertEquals(error.biasVarianceError, 0, biasVarianceTolerance);
        Assert.assertEquals(error.outlierProbabilityError, 0, outlierProbabilityTolerance);
    }

    @Test
    public void testResponsibility() {
        final SimpleInterval DUMMY = new SimpleInterval("dummy", 1, 2);
        final double VERY_LIKELY = 0.9;

        final AlleleFractionInitializer.Responsibility resp1 = new AlleleFractionInitializer.Responsibility(new AllelicCount(DUMMY, 1, 0));
        Assert.assertTrue(resp1.get(AlleleFractionIndicator.ALT_MINOR) > resp1.get(AlleleFractionIndicator.REF_MINOR));
        Assert.assertEquals(resp1.get(AlleleFractionIndicator.OUTLIER), AlleleFractionInitializer.INITIAL_OUTLIER_PROBABILITY);
        Assert.assertEquals(resp1.get(AlleleFractionIndicator.ALT_MINOR) + resp1.get(AlleleFractionIndicator.REF_MINOR) + resp1.get(AlleleFractionIndicator.OUTLIER), 1.0);

        final AlleleFractionInitializer.Responsibility resp2 = new AlleleFractionInitializer.Responsibility(new AllelicCount(DUMMY, 0, 1));
        Assert.assertTrue(resp2.get(AlleleFractionIndicator.ALT_MINOR) < resp2.get(AlleleFractionIndicator.REF_MINOR));
        Assert.assertEquals(resp1.get(AlleleFractionIndicator.OUTLIER), AlleleFractionInitializer.INITIAL_OUTLIER_PROBABILITY);

        final AlleleFractionInitializer.Responsibility resp3 = new AlleleFractionInitializer.Responsibility(new AllelicCount(DUMMY, 50, 150));
        Assert.assertTrue(resp3.get(AlleleFractionIndicator.REF_MINOR) > VERY_LIKELY);

        final AlleleFractionInitializer.Responsibility resp4 = new AlleleFractionInitializer.Responsibility(new AllelicCount(DUMMY, 150, 50));
        Assert.assertTrue(resp4.get(AlleleFractionIndicator.ALT_MINOR) > VERY_LIKELY);

        final AlleleFractionInitializer.Responsibility resp5 = new AlleleFractionInitializer.Responsibility(new AllelicCount(DUMMY, 1, 0));

        resp5.setFromLogPosteriors(1, 1, 1);
        Assert.assertEquals(resp5.get(AlleleFractionIndicator.ALT_MINOR), resp5.get(AlleleFractionIndicator.REF_MINOR));
        Assert.assertEquals(resp1.get(AlleleFractionIndicator.ALT_MINOR) + resp1.get(AlleleFractionIndicator.REF_MINOR) + resp1.get(AlleleFractionIndicator.OUTLIER), 1.0);

        resp5.setFromLogPosteriors(Double.NEGATIVE_INFINITY, 1, Double.NEGATIVE_INFINITY);
        Assert.assertEquals(resp5.get(AlleleFractionIndicator.REF_MINOR), 1.0);
        Assert.assertEquals(resp1.get(AlleleFractionIndicator.ALT_MINOR) + resp1.get(AlleleFractionIndicator.REF_MINOR) + resp1.get(AlleleFractionIndicator.OUTLIER), 1.0);

        resp5.setFromLogPosteriors(Math.log(1.0), Math.log(3.0), Double.NEGATIVE_INFINITY);
        Assert.assertEquals(resp5.get(AlleleFractionIndicator.ALT_MINOR), 1.0/4.0);
        Assert.assertEquals(resp1.get(AlleleFractionIndicator.ALT_MINOR) + resp1.get(AlleleFractionIndicator.REF_MINOR) + resp1.get(AlleleFractionIndicator.OUTLIER), 1.0);
    }
}