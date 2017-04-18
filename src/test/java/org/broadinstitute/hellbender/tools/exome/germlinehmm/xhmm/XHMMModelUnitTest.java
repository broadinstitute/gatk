package org.broadinstitute.hellbender.tools.exome.germlinehmm.xhmm;

import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriState;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriStateTransitionProbabilityCache;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.function.Function;

/**
 * Unit tests for {@link XHMMModel}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class XHMMModelUnitTest {

    private static double[] TEST_COVERAGE_VALUES = new double[] {
            100.0, 10.0, 3.0, 1.0, 0.5, 0.1, 0.00001, 0.0, -0.00001, -0.1, -0.5, -1.0, -3.0, -10.0, -100.0 };

    private static int[] TEST_TARGET_DISTANCES = new int[] {
            0, 10, 30, 100, 300, 1000, 3000, 10000, 300000, 1000000, 10000000};

    private static final double EPSILON = 1e-10;

    @Test(dataProvider = "testData")
    public void testInstantiation(final double eventStartProbability, final double meanEventSize,
                                  final double deletionMean, final double duplicationMean) {
         new XHMMModel(eventStartProbability, meanEventSize, deletionMean, duplicationMean);
    }

    @Test(dataProvider = "testData", dependsOnMethods = "testInstantiation")
    public void testHiddenStates(final double eventStartProbability, final double meanEventSize,
                               final double deletionMean, final double duplicationMean) {
        final XHMMModel model = new XHMMModel(eventStartProbability,
                meanEventSize, deletionMean, duplicationMean);

        final List<CopyNumberTriState> hiddenStates = model.hiddenStates();
        Assert.assertNotNull(hiddenStates);
        Assert.assertEquals(hiddenStates.size(), CopyNumberTriState.values().length);
        Assert.assertEquals(new HashSet<>(hiddenStates), new HashSet<>(Arrays.asList(CopyNumberTriState.values())));
    }

    @Test(dataProvider = "testData", dependsOnMethods = "testInstantiation")
    public void testLogEmissionProbability(final double eventStartProbability, final double meanEventSize,
                                           final double deletionMean, final double duplicationMean) {

        final XHMMModel model = new XHMMModel(eventStartProbability,
                meanEventSize, deletionMean, duplicationMean);
        final Target target = new Target("NAME");
        final double logDenominator = Math.log(Math.sqrt(2 * Math.PI));
        final Function<Double, Double> neutralEmission = x -> -.5 * (x * x) - logDenominator;
        final Function<Double, Double> deletionEmission = x -> -.5 * Math.pow(x - deletionMean, 2) - logDenominator;
        final Function<Double, Double> duplicationEmission = x -> -.5 * Math.pow(x - duplicationMean, 2) - logDenominator;
        for (final double coverage : TEST_COVERAGE_VALUES) {
            final double neutralObserved = model.logEmissionProbability(new XHMMEmissionData(coverage),
                    CopyNumberTriState.NEUTRAL, target);
            final double deletionObserved = model.logEmissionProbability(new XHMMEmissionData(coverage),
                    CopyNumberTriState.DELETION, target);
            final double duplicationObserved = model.logEmissionProbability(new XHMMEmissionData(coverage),
                    CopyNumberTriState.DUPLICATION, target);
            final double neutralExpected = neutralEmission.apply(coverage);
            final double deletionExpected = deletionEmission.apply(coverage);
            final double duplicationExpected = duplicationEmission.apply(coverage);
            Assert.assertEquals(neutralObserved, neutralExpected, "neutral emission for " + coverage);
            Assert.assertEquals(deletionObserved, deletionExpected);
            Assert.assertEquals(duplicationObserved, duplicationExpected);
        }
    }

    @Test(dataProvider = "testData", dependsOnMethods = "testInstantiation")
    public void testLogTransitionProbability(final double eventStartProbability, final double meanEventSize,
                                             final double deletionMean, final double duplicationMean) {
        final XHMMModel model = new XHMMModel(eventStartProbability,
                meanEventSize, deletionMean, duplicationMean);

        final Target fromTarget = new Target("target1", new SimpleInterval("1", 1000000000, 1000000100));
        for (final int distance : TEST_TARGET_DISTANCES) {
            final Target toTargetDownstream = new Target("target2", new SimpleInterval("1", 1000000100 + distance + 1, 1000000200 + distance + 1));
            final Target toTargetUpstream = new Target("target2", new SimpleInterval("1", 999999900 - distance - 1, 1000000000 - distance - 1));
            for (final Target toTarget : Arrays.asList(toTargetDownstream, toTargetUpstream)) {
                assertTransitionProbabilities(eventStartProbability, meanEventSize, model,  fromTarget, toTarget);
            }
        }
    }

    @Test(dataProvider = "testData", dependsOnMethods = "testInstantiation", expectedExceptions = IllegalArgumentException.class)
    public void testLogTransitionProbabilityWithNoDistance(final double eventStartProbability, final double meanEventSize,
                                                           final double deletionMean, final double duplicationMean) {
        final XHMMModel model = new XHMMModel(eventStartProbability,
                meanEventSize, deletionMean, duplicationMean);

        final Target fromTarget = new Target("target1");
        final Target toTarget = new Target("target2");

        assertTransitionProbabilities(eventStartProbability, meanEventSize, model, fromTarget, toTarget);
    }

    @Test(dataProvider = "testData", dependsOnMethods = "testInstantiation")
    public void testLogTransitionProbabilityOverlappingTargets(final double eventRate, final double meanEventSize,
                                                           final double deletionMean, final double duplicationMean) {

        final XHMMModel model = new XHMMModel(eventRate,
                meanEventSize, deletionMean, duplicationMean);
        final Target fromTarget = new Target("target1", new SimpleInterval("1", 1, 100));
        final Target toTarget = new Target("target2", new SimpleInterval("1", 50, 150));
        assertTransitionProbabilities(eventRate, meanEventSize, model, fromTarget, toTarget);
    }

    //the prior should be the infinite distance limit of the transition matrix
    @Test(dataProvider = "testData", dependsOnMethods = "testInstantiation")
    public void testLogPrior(final double eventStartProbability, final double meanEventSize,
                                             final double deletionMean, final double duplicationMean) {
        final XHMMModel model = new XHMMModel(eventStartProbability,
                meanEventSize, deletionMean, duplicationMean);
        final Target target = new Target("NAME");

        final CopyNumberTriStateTransitionProbabilityCache logTransitionProbabilityCache =
                new CopyNumberTriStateTransitionProbabilityCache(meanEventSize, eventStartProbability);

        for (final CopyNumberTriState state : CopyNumberTriState.values()) {
            Assert.assertEquals(model.logPriorProbability(state, target),
                    logTransitionProbabilityCache.logProbability(Integer.MAX_VALUE, state, CopyNumberTriState.NEUTRAL), 1e-10);
        }
    }

    // different chromosomes are independent Markov chains; hence the transition probabilities should equal the prior
    // probabilities
    @Test(dataProvider = "testData", dependsOnMethods = {"testInstantiation", "testLogPrior"})
    public void testLogTransitionProbabilityInDifferentChromosomes(final double eventStartProbability, final double meanEventSize,
                                                                   final double deletionMean, final double duplicationMean) {
        final XHMMModel model = new XHMMModel(eventStartProbability,
                meanEventSize, deletionMean, duplicationMean);

        final Target fromTarget = new Target("target1", new SimpleInterval("1", 1, 100));
        final Target toTarget = new Target("target2", new SimpleInterval("2", 1, 100));

        assertTransitionProbabilities(eventStartProbability, meanEventSize, model, fromTarget, toTarget);
    }

    @Test(dependsOnMethods = "testInstantiation")
    public void testDistance() {
        final XHMMModel model = new XHMMModel(0.5, 10, -1, 1);

        // different chromosomes
        Assert.assertEquals(XHMMModel.calculateDistance(new Target("target1", new SimpleInterval("1", 1, 100)),
                new Target("target2", new SimpleInterval("2", 1, 100))), Double.POSITIVE_INFINITY);

        // commutativity
        Assert.assertEquals(XHMMModel.calculateDistance(new Target("target1", new SimpleInterval("1", 1, 100)),
                new Target("target2", new SimpleInterval("1", 200, 300))),
                XHMMModel.calculateDistance(new Target("target1", new SimpleInterval("1", 200, 300)),
                        new Target("target2", new SimpleInterval("1", 1, 100))));

        Assert.assertEquals(XHMMModel.calculateDistance(new Target("target1", new SimpleInterval("1", 100, 200)),
                new Target("target2", new SimpleInterval("1", 200, 300))), 100, EPSILON);

        Assert.assertEquals(XHMMModel.calculateDistance(new Target("target1", new SimpleInterval("1", 100, 200)),
                new Target("target2", new SimpleInterval("1", 100, 110))), 45, EPSILON);

        Assert.assertEquals(XHMMModel.calculateDistance(new Target("target1", new SimpleInterval("1", 100, 200)),
                new Target("target2", new SimpleInterval("1", 150, 250))), 50, EPSILON);

    }


    @DataProvider(name = "testData")
    public Object[][] testData() {
        return new Object[][] {
                {1.0e-10D, 70000D, -1D, 1D},
                {1.0e-5D, 10000D, -3.3D, +3.3D},
                {.25D, 10D, -.5D, +.5D}
        };
    }

    // CopyNumberTriStateTransitionProbabilityCache is thoroughly tested.
    // Here we simply test that it is invoked  properly
    private void assertTransitionProbabilities(double eventStartProbability, double meanEventSize,
                                               XHMMModel model, Target fromTarget, Target toTarget) {

        final CopyNumberTriStateTransitionProbabilityCache logTransitionProbabilityCache =
                new CopyNumberTriStateTransitionProbabilityCache(meanEventSize, eventStartProbability);

        final double distance = XHMMModel.calculateDistance(fromTarget, toTarget);

        for (final CopyNumberTriState from : CopyNumberTriState.values()) {
            for (final CopyNumberTriState to : CopyNumberTriState.values()) {
                Assert.assertEquals(model.logTransitionProbability(from, fromTarget, to, toTarget),
                        logTransitionProbabilityCache.logProbability((int) distance, to, from), EPSILON);
            }
        }
    }

}
