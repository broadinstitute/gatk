package org.broadinstitute.hellbender.utils.recalibration;


import org.broadinstitute.hellbender.transformers.BQSRReadTransformer;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;


public final class RecalDatumUnitTest extends BaseTest {

    // --------------------------------------------------------------------------------
    //
    // merge case Provider
    //
    // --------------------------------------------------------------------------------

    private class RecalDatumTestProvider extends TestDataProvider {
        int exError, exTotal, reportedQual;

        private RecalDatumTestProvider(int E, int N, int reportedQual) {
            super(RecalDatumTestProvider.class);

            this.exError = E;
            this.exTotal = N;
            this.reportedQual = reportedQual;
        }

        public double getErrorRate() {
            return (exError + 1) / (1.0 * (exTotal + 2));
        }

        public int getReportedQual() {
            return reportedQual;
        }

        public RecalDatum makeRecalDatum() {
            return new RecalDatum((long)exTotal, (double)exError, (byte)getReportedQual());
        }

        @Override
        public String toString() {
            return String.format("exError=%d, exTotal=%d, reportedQual=%d", exError, exTotal, reportedQual);
        }
    }

    private static boolean createdDatumTestProviders = false;

    @DataProvider(name = "RecalDatumTestProvider")
    public Object[][] makeRecalDatumTestProvider() {
        if ( !createdDatumTestProviders ) {
            for ( int E : Arrays.asList(1, 10, 100, 1000, 10000) )
                for ( int N : Arrays.asList(10, 100, 1000, 10000, 100000, 1000000) )
                    for ( int reportedQual : Arrays.asList(10, 20) )
                        if ( E <= N )
                            new RecalDatumTestProvider(E, N, reportedQual);
            createdDatumTestProviders = true;
        }

        return RecalDatumTestProvider.getTests(RecalDatumTestProvider.class);
    }

    @Test(dataProvider = "RecalDatumTestProvider")
    public void testRecalDatumBasics(RecalDatumTestProvider cfg) {
        final RecalDatum datum = cfg.makeRecalDatum();
        assertBasicFeaturesOfRecalDatum(datum, cfg);
    }

    private static void assertBasicFeaturesOfRecalDatum(final RecalDatum datum, final RecalDatumTestProvider cfg) {
        Assert.assertEquals(datum.getNumMismatches(), cfg.exError, 1E-6);
        Assert.assertEquals(datum.getNumObservations(), cfg.exTotal, 1E-6);
        if ( cfg.getReportedQual() != -1 )
            Assert.assertEquals(datum.getEstimatedQReportedAsByte(), cfg.getReportedQual());
        assertEqualsDoubleSmart(datum.getEmpiricalErrorRate(), cfg.getErrorRate());

        final double e = datum.getEmpiricalQuality();
        Assert.assertTrue(datum.getEmpiricalQualityAsByte() >= Math.floor(e));
        Assert.assertTrue(datum.getEmpiricalQualityAsByte() <= Math.ceil(e));
        Assert.assertNotNull(datum.toString());
    }

    @Test(dataProvider = "RecalDatumTestProvider")
    public void testRecalDatumCopyAndCombine(RecalDatumTestProvider cfg) {
        final RecalDatum datum = cfg.makeRecalDatum();
        final RecalDatum copy = new RecalDatum(datum);
        assertBasicFeaturesOfRecalDatum(copy, cfg);

        RecalDatumTestProvider combinedCfg = new RecalDatumTestProvider(cfg.exError * 2, cfg.exTotal * 2, cfg.reportedQual);
        copy.combine(datum);
        assertBasicFeaturesOfRecalDatum(copy, combinedCfg);
    }

    @Test(dataProvider = "RecalDatumTestProvider")
    public void testRecalDatumModification(RecalDatumTestProvider cfg) {
        RecalDatum datum = cfg.makeRecalDatum();
        datum.setEmpiricalQuality(10.1);
        Assert.assertEquals(datum.getEmpiricalQuality(), 10.1);

        datum.setEstimatedQReported(10.1);
        Assert.assertEquals(datum.getEstimatedQReported(), 10.1);
        Assert.assertEquals(datum.getEstimatedQReportedAsByte(), 10);

        datum = cfg.makeRecalDatum();
        cfg.exTotal = 100000;
        datum.setNumObservations(cfg.exTotal);
        assertBasicFeaturesOfRecalDatum(datum, cfg);

        datum = cfg.makeRecalDatum();
        cfg.exError = 1000;
        datum.setNumMismatches(cfg.exError);
        assertBasicFeaturesOfRecalDatum(datum, cfg);

        datum = cfg.makeRecalDatum();
        datum.increment(true);
        cfg.exError++;
        cfg.exTotal++;
        assertBasicFeaturesOfRecalDatum(datum, cfg);

        datum = cfg.makeRecalDatum();
        datum.increment(false);
        cfg.exTotal++;
        assertBasicFeaturesOfRecalDatum(datum, cfg);

        datum = cfg.makeRecalDatum();
        datum.incrementNumObservations(2);
        cfg.exTotal += 2;
        assertBasicFeaturesOfRecalDatum(datum, cfg);

        datum = cfg.makeRecalDatum();
        datum.incrementNumMismatches(2);
        cfg.exError += 2;
        assertBasicFeaturesOfRecalDatum(datum, cfg);


        datum = cfg.makeRecalDatum();
        datum.increment(10, 5);
        cfg.exError += 5;
        cfg.exTotal += 10;
        assertBasicFeaturesOfRecalDatum(datum, cfg);
    }

    @Test
    public void testNoObs() {
        final RecalDatum rd = new RecalDatum(0L, 0.0, (byte)10);
        Assert.assertEquals(rd.getEmpiricalErrorRate(), 0.0);
    }

    @Test
    public void testlogQempPrior() {
        for ( int Qemp = 0; Qemp <= QualityUtils.MAX_SAM_QUAL_SCORE; Qemp++ ) {
            for ( int Qrep = 0; Qrep <= QualityUtils.MAX_SAM_QUAL_SCORE; Qrep++ ) {
                final double logPrior = RecalDatum.logQempPrior(Qemp, Qrep);
                Assert.assertTrue(logPrior < 0.0);
                Assert.assertFalse(Double.isInfinite(logPrior));
                Assert.assertFalse(Double.isNaN(logPrior));
            }
        }

        final int Qrep = 20;
        int maxQemp = -1;
        double maxQempValue = -Double.MAX_VALUE;
        for ( int Qemp = 0; Qemp <= QualityUtils.MAX_SAM_QUAL_SCORE; Qemp++ ) {
            final double logPrior = RecalDatum.logQempPrior(Qemp, Qrep);
            if ( logPrior > maxQempValue ) {
                maxQemp = Qemp;
                maxQempValue = logPrior;
            }
        }
        Assert.assertEquals(maxQemp, Qrep);
    }

    @Test
    public void testBayesianEstimateOfEmpiricalQuality() {

        final int Qrep = 20;

        // test no shift
        Assert.assertEquals(RecalDatum.bayesianEstimateOfEmpiricalQuality(0, 0, Qrep), (double) Qrep);
        Assert.assertEquals(RecalDatum.bayesianEstimateOfEmpiricalQuality(10, 0, Qrep), (double) Qrep);
        Assert.assertEquals(RecalDatum.bayesianEstimateOfEmpiricalQuality(1000, 10, Qrep), (double) Qrep);

        // test small shift
        Assert.assertEquals(RecalDatum.bayesianEstimateOfEmpiricalQuality(10, 10, Qrep), Qrep - 1.0);
        Assert.assertEquals(RecalDatum.bayesianEstimateOfEmpiricalQuality(1000, 0, Qrep), Qrep + 1.0);

        // test medium shift
        Assert.assertEquals(RecalDatum.bayesianEstimateOfEmpiricalQuality(10000, 0, Qrep), Qrep + 3.0);
        Assert.assertEquals(RecalDatum.bayesianEstimateOfEmpiricalQuality(10000, 10, Qrep), Qrep + 3.0);

        // test large shift
        Assert.assertEquals(RecalDatum.bayesianEstimateOfEmpiricalQuality(100000, 10, Qrep), Qrep + 8.0);
        Assert.assertEquals(RecalDatum.bayesianEstimateOfEmpiricalQuality(1000000, 10, Qrep), Qrep + 16.0);
    }

    @Test
    public void testlogQempLikelihood() {

        final double[] Qemps = new double[] { 0.0, 10.0, 20.0, 30.0 };
        final int[] observations = new int[] {0, 10, 1000, 1000000};
        final int[] errors = new int[] {0, 10, 1000, 1000000};

        for ( double Qemp : Qemps ) {
            for ( int observation : observations ) {
                for ( int error : errors ) {
                    if ( error > observation )
                        continue;

                    final double loglikelihood = RecalDatum.logQempLikelihood(Qemp, observation, error);
                    Assert.assertTrue(observation == 0 ? MathUtils.compareDoubles(loglikelihood, 0.0) == 0 : loglikelihood < 0.0);
                    Assert.assertFalse(Double.isInfinite(loglikelihood));
                    Assert.assertFalse(Double.isNaN(loglikelihood));
                }
            }
        }

        long bigNum = new Long((long) Integer.MAX_VALUE);
        bigNum *= 2L;
        final double loglikelihood = RecalDatum.logQempLikelihood(30, bigNum, 100000);
        Assert.assertTrue(loglikelihood < 0.0);
        Assert.assertFalse(Double.isInfinite(loglikelihood));
        Assert.assertFalse(Double.isNaN(loglikelihood));
    }

    @Test
    public void basicHierarchicalBayesianQualityEstimateTest() {

        for( double epsilon = 15.0; epsilon <= 60.0; epsilon += 2.0 ) {
            double RG_Q = 45.0;
            RecalDatum RG = new RecalDatum( (long)100000000, (long) (100000000 * 1.0 / (Math.pow(10.0, RG_Q / 10.0))), (byte)RG_Q);
            double Q = 30.0;
            RecalDatum QS = new RecalDatum( (long)100000000, (long) (100000000 * 1.0 / (Math.pow(10.0, Q / 10.0))), (byte)Q);
            RecalDatum COV = new RecalDatum( (long)15, (long) 1, (byte)45.0); // no data here so Bayesian prior has a huge effect on the empirical quality

            // initial epsilon condition shouldn't matter when there are a lot of observations
            Assert.assertEquals(BQSRReadTransformer.hierarchicalBayesianQualityEstimate(epsilon, RG, QS, COV), Q, 1E-4);
        }

        for( double epsilon = 15.0; epsilon <= 60.0; epsilon += 2.0 ) {
            double RG_Q = 45.0;
            RecalDatum RG = new RecalDatum( (long)10, (long) (10 * 1.0 / (Math.pow(10.0, RG_Q / 10.0))), (byte)RG_Q);
            double Q = 30.0;
            RecalDatum QS = new RecalDatum( (long)10, (long) (10 * 1.0 / (Math.pow(10.0, Q / 10.0))), (byte)Q);
            RecalDatum COV = new RecalDatum( (long)15, (long) 1, (byte)45.0); // no data here so Bayesian prior has a huge effect on the empirical quality

            // initial epsilon condition dominates when there is no data
            Assert.assertEquals(BQSRReadTransformer.hierarchicalBayesianQualityEstimate(epsilon, RG, QS, COV), epsilon, 1E-4);
        }

    }
}