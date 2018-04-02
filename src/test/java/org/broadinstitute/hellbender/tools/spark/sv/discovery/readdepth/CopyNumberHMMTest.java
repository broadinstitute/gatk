package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.RealVector;
import org.broadinstitute.hellbender.utils.hmm.ViterbiAlgorithm;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

public class CopyNumberHMMTest extends BaseTest {

    private static final int NUM_STATES = 5;
    private static final double TRANS_PROB = 0.01;

    @Test(groups = "sv")
    public void testUniformPrior() {
        final RealVector prior = CopyNumberHMM.uniformPrior(NUM_STATES);
        Assert.assertEquals(prior.getDimension(), NUM_STATES);
        for (final double val : prior.toArray()) {
            Assert.assertEquals(val, 1.0 / NUM_STATES);
        }
    }

    @Test(groups = "sv")
    public void testPositionsList() {
        final List<Integer> positionsList = CopyNumberHMM.positionsList(10);
        for (int i = 0; i < 10; i++) {
            Assert.assertEquals(positionsList.get(i).intValue(), i);
        }
    }

    @Test(groups = "sv")
    public void testHiddenStates() {
        final RealVector prior = CopyNumberHMM.uniformPrior(NUM_STATES);
        final CopyNumberHMM hmm = new CopyNumberHMM(prior, TRANS_PROB);
        final List<Integer> hiddenStates = hmm.hiddenStates();
        Assert.assertEquals(hiddenStates.size(), NUM_STATES);
    }

    @Test(groups = "sv")
    public void testLogPriorProbability() {
        final RealVector prior = CopyNumberHMM.uniformPrior(NUM_STATES);
        final CopyNumberHMM hmm = new CopyNumberHMM(prior, TRANS_PROB);
        for (int i = 0; i < NUM_STATES; i++) {
            Assert.assertEquals(hmm.logPriorProbability(i, 0), Math.log(prior.getEntry(i)));
        }
    }

    @Test(groups = "sv")
    public void testLogTransitionProbability() {
        final RealVector prior = CopyNumberHMM.uniformPrior(NUM_STATES);
        final CopyNumberHMM hmm = new CopyNumberHMM(prior, TRANS_PROB);
        final double logTransProb = Math.log(TRANS_PROB / (NUM_STATES - 1));
        final double logNoTransProb = Math.log(1.0 - TRANS_PROB);
        for (int i = 0; i < NUM_STATES; i++) {
            for (int j = 0; j < NUM_STATES; j++) {
                final double expectedTransProb = (i == j) ? logNoTransProb : logTransProb;
                Assert.assertEquals(hmm.logTransitionProbability(i, 0, j, 1), expectedTransProb);
            }
        }
    }

    @Test(groups = "sv")
    public void testLogEmissionProbability() {
        final RealVector prior = CopyNumberHMM.uniformPrior(NUM_STATES);
        final CopyNumberHMM hmm = new CopyNumberHMM(prior, TRANS_PROB);
        for (int i = 0; i < NUM_STATES; i++) {
            final NormalDistribution emissionDist = new NormalDistribution(0.5 * i, CopyNumberHMM.COPY_RATIO_STDEV);
            for (int j = 0; j < 5; j++) {
                final double copyRatio = 0.5 * j;
                Assert.assertEquals(hmm.logEmissionProbability(copyRatio, i, 0), emissionDist.logDensity(copyRatio));
            }
        }
    }

    @DataProvider(name = "viterbiTestData")
    public Object[][] getViterbiData() {
        return new Object[][]{
                {Arrays.asList(1.1, 1.02, 0.95, 1.03, 1.0, 0.87, 1.22, 1.05, 1.04, 1.37, 0.9, 1.15), Arrays.asList(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)},
                {Arrays.asList(1.4, 1.5, 1.8, 1.6, 1.5, 1.4, 1.45, 1.65, 1.3, 1.6, 1.5), Arrays.asList(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)},
                {Arrays.asList(0.1, 0.2, 0.1, 0.05, 0.0, 0.3, 0.2, 0.05, 0.0, 0.0, 0.25), Arrays.asList(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)},
                {Arrays.asList(1.1, 0.9, 0.8, 0.9, 1.2, 1.5, 1.6, 1.4, 1.5, 1.5, 1.6), Arrays.asList(2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3)}
        };
    }
    @Test(groups = "sv",
            dataProvider = "viterbiTestData")
    public void testViterbi(final List<Double> copyRatios, final List<Integer> expectedStates) {
        final RealVector prior = CopyNumberHMM.uniformPrior(NUM_STATES);
        final CopyNumberHMM hmm = new CopyNumberHMM(prior, TRANS_PROB);

        final List<Integer> positions = CopyNumberHMM.positionsList(copyRatios.size());
        final List<Integer> states = ViterbiAlgorithm.apply(copyRatios, positions, hmm);
        Assert.assertEquals(states, expectedStates);
    }
}