package org.broadinstitute.hellbender.tools.exome.hmm;

import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.function.Function;

/**
 * Unit tests for {@link CopyNumberTriStateHiddenMarkovModel}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class CopyNumberTriStateHiddenMarkovModelUnitTest {

    private static double[] TEST_COVERAGE_VALUES = new double[] {
            100.0, 10.0, 3.0, 1.0, 0.5, 0.1, 0.00001, 0.0, -0.00001, -0.1, -0.5, -1.0, -3.0, -10.0, -100.0 };

    private static int[] TEST_TARGET_DISTANCES = new int[] {
            0, 10, 30, 100, 300, 1000, 3000, 10000, 300000, 1000000, 10000000};

    @Test(dataProvider = "testData")
    public void testInstantiation(final double eventRate, final double meanNumberOfTargetsPerEvent,
                                  final double deletionMean, final double duplicationMean,
                                  final double averageDistanceBetweenTargetsInEvent) {
         new CopyNumberTriStateHiddenMarkovModel(eventRate,
                meanNumberOfTargetsPerEvent, deletionMean, duplicationMean, averageDistanceBetweenTargetsInEvent);
    }

    @Test(dataProvider = "testData", dependsOnMethods = "testInstantiation")
    public void testHiddenStates(final double eventRate, final double meanNumberOfTargetsPerEvent,
                               final double deletionMean, final double duplicationMean,
                               final double averageDistanceBetweenTargetsInEvent) {

        final CopyNumberTriStateHiddenMarkovModel model = new CopyNumberTriStateHiddenMarkovModel(eventRate,
                meanNumberOfTargetsPerEvent, deletionMean, duplicationMean, averageDistanceBetweenTargetsInEvent);

        final List<CopyNumberTriState> hiddenStates = model.hiddenStates();
        Assert.assertNotNull(hiddenStates);
        Assert.assertEquals(hiddenStates.size(), CopyNumberTriState.values().length);
        Assert.assertEquals(new HashSet<>(hiddenStates), new HashSet<>(Arrays.asList(CopyNumberTriState.values())));
    }

    @Test(dataProvider = "testData", dependsOnMethods = "testInstantiation")
    public void testLogEmissionProbability(final double eventRate, final double meanNumberOfTargetsPerEvent,
                                           final double deletionMean, final double duplicationMean,
                                           final double averageDistanceBetweenTargetsInEvent) {

        final CopyNumberTriStateHiddenMarkovModel model = new CopyNumberTriStateHiddenMarkovModel(eventRate,
                meanNumberOfTargetsPerEvent, deletionMean, duplicationMean, averageDistanceBetweenTargetsInEvent);
        final Target target = new Target("NAME");
        final double logDenominator = Math.log(Math.sqrt(2 * Math.PI));
        final Function<Double, Double> neutralEmission = x -> -.5 * (x * x) - logDenominator;
        final Function<Double, Double> deletionEmission = x -> -.5 * Math.pow(x - deletionMean, 2) - logDenominator;
        final Function<Double, Double> duplicationEmission = x -> -.5 * Math.pow(x - duplicationMean, 2) - logDenominator;
        for (final double coverage : TEST_COVERAGE_VALUES) {
            final double neutralObserved = model.logEmissionProbability(coverage, CopyNumberTriState.NEUTRAL, target);
            final double deletionObserved = model.logEmissionProbability(coverage, CopyNumberTriState.DELETION, target);
            final double duplicationObserved = model.logEmissionProbability(coverage, CopyNumberTriState.DUPLICATION, target);
            final double neutralExpected = neutralEmission.apply(coverage);
            final double deletionExpected = deletionEmission.apply(coverage);
            final double duplicationExpected = duplicationEmission.apply(coverage);
            Assert.assertEquals(neutralObserved, neutralExpected, "neutral emission for " + coverage);
            Assert.assertEquals(deletionObserved, deletionExpected);
            Assert.assertEquals(duplicationObserved, duplicationExpected);
        }
    }

    @Test(dataProvider = "testData", dependsOnMethods = "testInstantiation")
    public void testLogTransitionProbability(final double eventRate, final double meanNumberOfTargetsPerEvent,
                                             final double deletionMean, final double duplicationMean,
                                             final double averageDistanceBetweenTargetsInEvent) {

        final CopyNumberTriStateHiddenMarkovModel model = new CopyNumberTriStateHiddenMarkovModel(eventRate,
                meanNumberOfTargetsPerEvent, deletionMean, duplicationMean, averageDistanceBetweenTargetsInEvent);
        final Target fromTarget = new Target("target1", new SimpleInterval("1", 1000000000, 1000000100));
        for (final int distance : TEST_TARGET_DISTANCES) {
            final Target toTargetDownstream = new Target("target2", new SimpleInterval("1", 1000000100 + distance + 1, 1000000200 + distance + 1));
            final Target toTargetUpstream = new Target("target2", new SimpleInterval("1", 999999900 - distance - 1, 1000000000 - distance - 1));
            for (final Target toTarget : Arrays.asList(toTargetDownstream, toTargetUpstream)) {
                assertTransitionProbabilities(eventRate, meanNumberOfTargetsPerEvent, averageDistanceBetweenTargetsInEvent, model, distance, fromTarget, toTarget);
            }
        }
    }

    @Test(dataProvider = "testData", dependsOnMethods = "testInstantiation")
    public void testLogTransitionProbabilityWithNoDistance(final double eventRate, final double meanNumberOfTargetsPerEvent,
                                                           final double deletionMean, final double duplicationMean,
                                                           final double averageDistanceBetweenTargetsInEvent) {

        final CopyNumberTriStateHiddenMarkovModel model = new CopyNumberTriStateHiddenMarkovModel(eventRate,
                meanNumberOfTargetsPerEvent, deletionMean, duplicationMean, averageDistanceBetweenTargetsInEvent);
        final double distance = 0.0;
        final Target fromTarget = new Target("target1");
        final Target toTarget = new Target("target2");
        assertTransitionProbabilities(eventRate, meanNumberOfTargetsPerEvent, averageDistanceBetweenTargetsInEvent, model, distance, fromTarget, toTarget);
    }

    @Test(dataProvider = "testData", dependsOnMethods = "testInstantiation")
    public void testLogTransitionProbabilityInDifferentChromosomes(final double eventRate, final double meanNumberOfTargetsPerEvent,
                                                           final double deletionMean, final double duplicationMean,
                                                           final double averageDistanceBetweenTargetsInEvent) {

        final CopyNumberTriStateHiddenMarkovModel model = new CopyNumberTriStateHiddenMarkovModel(eventRate,
                meanNumberOfTargetsPerEvent, deletionMean, duplicationMean, averageDistanceBetweenTargetsInEvent);
        final double distance = Double.POSITIVE_INFINITY;
        final Target fromTarget = new Target("target1", new SimpleInterval("1", 1, 100));
        final Target toTarget = new Target("target2", new SimpleInterval("2", 1, 100));
        assertTransitionProbabilities(eventRate, meanNumberOfTargetsPerEvent, averageDistanceBetweenTargetsInEvent, model, distance, fromTarget, toTarget);
    }

    private void assertTransitionProbabilities(double eventRate, double meanNumberOfTargetsPerEvent, double averageDistanceBetweenTargetsInEvent, CopyNumberTriStateHiddenMarkovModel model, double distance, Target fromTarget, Target toTarget) {
        final double logP = Math.log(eventRate);
        final double log1Minus2P = Math.log1p(-2 * eventRate);
        final double logF = -distance / averageDistanceBetweenTargetsInEvent;
        final double logQ = -Math.log(meanNumberOfTargetsPerEvent);
        final double log1MinusF = Math.log1p(-Math.exp(logF));
        final double observedDistance = model.calculateDistance(fromTarget, toTarget);
        final double log1MinusQ = Math.log1p(-(1.0 / meanNumberOfTargetsPerEvent));
        //logF + log1MinusQ, log1MinusF + logP

        Assert.assertEquals(observedDistance, distance);
        Assert.assertEquals(model.logTransitionProbability(CopyNumberTriState.NEUTRAL, fromTarget,
                CopyNumberTriState.DELETION, toTarget), logP, " " + eventRate);
        Assert.assertEquals(model.logTransitionProbability(CopyNumberTriState.NEUTRAL, fromTarget,
                CopyNumberTriState.DUPLICATION, toTarget), logP);
        Assert.assertEquals(model.logTransitionProbability(CopyNumberTriState.NEUTRAL, fromTarget,
                CopyNumberTriState.NEUTRAL, toTarget), log1Minus2P, "eventRate = " + eventRate);
        Assert.assertEquals(model.logTransitionProbability(CopyNumberTriState.DELETION, fromTarget,
                        CopyNumberTriState.DUPLICATION, toTarget), log1MinusF + logP,
                "avgDist = " + averageDistanceBetweenTargetsInEvent + " distance = " + distance + " eventRate = " + eventRate);
        Assert.assertEquals(model.logTransitionProbability(CopyNumberTriState.DUPLICATION, fromTarget,
                        CopyNumberTriState.DELETION, toTarget), log1MinusF + logP,
                "avgDist = " + averageDistanceBetweenTargetsInEvent + " distance = " + distance + " eventRate = " + eventRate);
        Assert.assertEquals(model.logTransitionProbability(CopyNumberTriState.DELETION, fromTarget,
                        CopyNumberTriState.NEUTRAL, toTarget), GATKProtectedMathUtils.naturalLogSumExp(logF + logQ, log1MinusF + log1Minus2P),
                epsilon(model.logTransitionProbability(CopyNumberTriState.DELETION, fromTarget,
                        CopyNumberTriState.NEUTRAL, toTarget), GATKProtectedMathUtils.naturalLogSumExp(logF + logQ, log1MinusF + log1Minus2P)));
        Assert.assertEquals(model.logTransitionProbability(CopyNumberTriState.DUPLICATION, fromTarget,
                        CopyNumberTriState.NEUTRAL, toTarget), GATKProtectedMathUtils.naturalLogSumExp(logF + logQ, log1MinusF + log1Minus2P),
                epsilon(model.logTransitionProbability(CopyNumberTriState.DUPLICATION, fromTarget,
                        CopyNumberTriState.NEUTRAL, toTarget), GATKProtectedMathUtils.naturalLogSumExp(logF + logQ, log1MinusF + log1Minus2P)));
        Assert.assertEquals(model.logTransitionProbability(CopyNumberTriState.DUPLICATION, fromTarget,
                        CopyNumberTriState.DUPLICATION, toTarget), GATKProtectedMathUtils.naturalLogSumExp(logF + log1MinusQ, log1MinusF + logP),
                epsilon(model.logTransitionProbability(CopyNumberTriState.DUPLICATION, fromTarget,
                        CopyNumberTriState.DUPLICATION, toTarget), GATKProtectedMathUtils.naturalLogSumExp(logF + log1MinusQ, log1MinusF + logP)),
                "T = " + meanNumberOfTargetsPerEvent + " avgDist = " + averageDistanceBetweenTargetsInEvent + " distance = " + distance + " eventRate = " + eventRate);
        Assert.assertEquals(model.logTransitionProbability(CopyNumberTriState.DELETION, fromTarget,
                        CopyNumberTriState.DELETION, toTarget), GATKProtectedMathUtils.naturalLogSumExp(logF + log1MinusQ, log1MinusF + logP),
                epsilon(model.logTransitionProbability(CopyNumberTriState.DELETION, fromTarget,
                        CopyNumberTriState.DELETION, toTarget), GATKProtectedMathUtils.naturalLogSumExp(    logF + log1MinusQ, log1MinusF + logP)),
                "T = " + meanNumberOfTargetsPerEvent + " avgDist = " + averageDistanceBetweenTargetsInEvent + " distance = " + distance + " eventRate = " + eventRate);
    }

    @Test(dataProvider = "testData", dependsOnMethods = "testInstantiation")
    public void testLogTransitionProbabilityOverlappingTargets(final double eventRate, final double meanNumberOfTargetsPerEvent,
                                                           final double deletionMean, final double duplicationMean,
                                                           final double averageDistanceBetweenTargetsInEvent) {

        final CopyNumberTriStateHiddenMarkovModel model = new CopyNumberTriStateHiddenMarkovModel(eventRate,
                meanNumberOfTargetsPerEvent, deletionMean, duplicationMean, averageDistanceBetweenTargetsInEvent);
        final double distance = 0.0;
        final Target fromTarget = new Target("target1", new SimpleInterval("1", 1, 100));
        final Target toTarget = new Target("target2", new SimpleInterval("1", 50, 150));
        assertTransitionProbabilities(eventRate, meanNumberOfTargetsPerEvent, averageDistanceBetweenTargetsInEvent, model, distance, fromTarget, toTarget);
    }

    @Test(dataProvider = "testData", dependsOnMethods = "testInstantiation")
    public void testLogPrior(final double eventRate, final double meanNumberOfTargetsPerEvent,
                                             final double deletionMean, final double duplicationMean,
                                             final double averageDistanceBetweenTargetsInEvent) {
        final CopyNumberTriStateHiddenMarkovModel model = new CopyNumberTriStateHiddenMarkovModel(eventRate,
                meanNumberOfTargetsPerEvent, deletionMean, duplicationMean, averageDistanceBetweenTargetsInEvent);
        final Target target = new Target("NAME");
        Assert.assertEquals(model.logPriorProbability(CopyNumberTriState.DELETION, target), Math.log(eventRate));
        Assert.assertEquals(model.logPriorProbability(CopyNumberTriState.DUPLICATION, target), Math.log(eventRate));
        Assert.assertEquals(model.logPriorProbability(CopyNumberTriState.NEUTRAL, target), Math.log1p(- 2*eventRate));
    }

    private static double epsilon(final double a, final double b) {
        return 0.000001 * Math.min(Math.abs(a), Math.abs(b));
    }
    @DataProvider(name = "testData")
    public Object[][] testData() {
        return new Object[][] {
                {1.0e-10D, 6.0D, -1D, 1D, 70000D},
                {1.0e-5D,  3.0D, -3.3D, +3.3D, 10000D},
                {.25D,      100.0D, -.5D, +.5D, 10D}
        };
    }

}
