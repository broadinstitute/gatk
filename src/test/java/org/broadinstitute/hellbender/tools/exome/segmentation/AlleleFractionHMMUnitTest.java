package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionGlobalParameters;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.hmm.ForwardBackwardAlgorithm;
import org.broadinstitute.hellbender.utils.hmm.ViterbiAlgorithm;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionHMMUnitTest {
    private static final AlleleFractionGlobalParameters NO_BIAS_OR_OUTLIERS_PARAMS = new AlleleFractionGlobalParameters(1.0, 1e-10, 1e-10);

    @Test
    public void constructorTest() {
        final double memoryLength = 1e6;
        final List<Double> minorAlleleFractions = Arrays.asList(0.1, 0.5, 0.23);
        final AlleleFractionGlobalParameters params = new AlleleFractionGlobalParameters(0.1, 0.01, 0.03);
        final AlleleFractionHMM model = new AlleleFractionHMM(minorAlleleFractions,
                memoryLength);

        Assert.assertEquals(memoryLength, model.getMemoryLength());
        for (int n = 0; n < minorAlleleFractions.size(); n++) {
            Assert.assertEquals(minorAlleleFractions.get(n), model.getMinorAlleleFraction(n));
        }
    }

    // if all states have the same minor fraction, then regardless of data the hidden state probabilities are equal
    @Test
    public void equalMinorFractionsTest() {
        final List<Double> minorAlleleFractions = Arrays.asList(0.3, 0.3, 0.3);
        final double memoryLength = 1e3;
        final AlleleFractionGlobalParameters params = new AlleleFractionGlobalParameters(0.1, 0.01, 0.03);
        final AlleleFractionHMM model = new AlleleFractionHMM(minorAlleleFractions,
                memoryLength);

        final Random random = new Random(13);
        final int chainLength = 10000;
        final List<SimpleInterval> positions = new ArrayList<>();
        final List<AllelicCount> data = new ArrayList<>();
        int position = 1;
        for (int n = 0; n < chainLength; n++) {
            position += random.nextInt((int) (2*memoryLength));
            final SimpleInterval interval = new SimpleInterval("chr1", position, position);
            positions.add(interval);
            data.add(new AllelicCount(interval, random.nextInt(30) + 1, random.nextInt(30) + 1));
        }

        final ForwardBackwardAlgorithm.Result<AllelicCount, SimpleInterval, Integer> fbResult =
                ForwardBackwardAlgorithm.apply(data, positions, model);

        for (int pos = 0; pos < chainLength; pos++) {
            for (int state = 0; state < minorAlleleFractions.size(); state++) {
                Assert.assertEquals(fbResult.logProbability(pos, state), -Math.log(minorAlleleFractions.size()), 1e-5);
            }
        }
    }

    // if we ignore ref bias, choose  well-separated minor fractions e.g. 0.1 and 0.5, and give really obvious
    // allele counts like 10/90 (obviously 0.1) and 50/50 (obviously 0.5), the Viterbi algorithm should make the right calls
    // this essentially tests whether we correctly defined the likelihood in terms of methods in the Allele Fraction model
    @Test
    public void obviousCallsTest() {
        final List<Double> minorAlleleFractions = Arrays.asList(0.1, 0.5);
        final double memoryLength = 1e3;
        final AlleleFractionHMM model = new AlleleFractionHMM(minorAlleleFractions,
                memoryLength);

        final Random random = new Random(13);
        final int chainLength = 10000;
        final List<SimpleInterval> positions = new ArrayList<>();
        final List<AllelicCount> data = new ArrayList<>();
        int position = 1;
        for (int n = 0; n < chainLength; n++) {
            position += random.nextInt((int) (2*memoryLength)) + (int) memoryLength;
            final SimpleInterval interval = new SimpleInterval("chr1", position, position);
            positions.add(interval);

            if (n < 2500) {
                data.add(new AllelicCount(interval, 10, 90));   // minor fraction = 0.1
            } else if ( n < 5000) {
                data.add(new AllelicCount(interval, 50, 50));   // minor fraction = 0.5
            } else if ( n < 7500) {
                data.add(new AllelicCount(interval, 90, 10));   // minor fraction = 0.1
            } else {
                data.add(new AllelicCount(interval, 40, 60));   // minor fraction = 0.5
            }
        }

        final List<Integer> states = ViterbiAlgorithm.apply(data, positions, model);

        for (int pos = 0; pos < chainLength; pos++) {
            if (pos < 2500) {
                Assert.assertEquals((int) states.get(pos), 0);
            } else if (pos < 5000) {
                Assert.assertEquals((int) states.get(pos), 1);
            } else if (pos < 7500) {
                Assert.assertEquals((int) states.get(pos), 0);
            } else {
                Assert.assertEquals((int) states.get(pos), 1);
            }
        }
    }

    // if weights are equal, transition probabilities should equal prior probabilities far past the
    // memory length and should forbid transitions at zero distance
    @Test
    public void testTransitionProbabilities() {
        final List<Double> minorAlleleFractions = Arrays.asList(0.1, 0.2, 0.3, 0.4, 0.23, 0.11);
        final double memoryLength = 5e6;
        final AlleleFractionHMM model = new AlleleFractionHMM(minorAlleleFractions,
                memoryLength);

        final int pos1 = 100;
        final int pos2 = (int) (pos1 + 100*memoryLength);  // really far!!!
        final int pos3 = (int) (pos1 + memoryLength/1000);  // really close!!!
        final SimpleInterval position1 = new SimpleInterval("chr1", pos1, pos1);
        final SimpleInterval position2 = new SimpleInterval("chr1", pos2, pos2);
        final SimpleInterval position3 = new SimpleInterval("chr1", pos3, pos3);

        for (int fromState = 0; fromState < minorAlleleFractions.size(); fromState++) {
            for (int toState = 0; toState < minorAlleleFractions.size(); toState++) {
                // long distances
                Assert.assertEquals(model.logTransitionProbability(fromState, position1, toState, position2),
                        model.logPriorProbability(toState, position2), 1e-3);

                //short distances
                if (toState == fromState) {
                    Assert.assertEquals(model.logTransitionProbability(fromState, position1, toState, position3),
                            Math.log(1), 1e-3);
                } else {
                    Assert.assertTrue(model.logTransitionProbability(fromState, position1, toState, position3) < Math.log(1e-3));
                }
            }
        }
    }

    // test that the prior probabilities are flat
    @Test
    public void logPriorTest() {
        final List<Double> minorAlleleFractions = Arrays.asList(0.1, 0.2, 0.3, 0.4, 0.5, 0.33);
        final double memoryLength = 5e6;
        final AlleleFractionHMM model = new AlleleFractionHMM(minorAlleleFractions,
                memoryLength);

        final SimpleInterval position1 = new SimpleInterval("chr1", 100, 100);
        final SimpleInterval position2 = new SimpleInterval("chr2", 20000, 20000);
        for (int n = 0; n < minorAlleleFractions.size(); n++) {
            Assert.assertEquals(model.logPriorProbability(n, position1), -Math.log(minorAlleleFractions.size()));
            Assert.assertEquals(model.logPriorProbability(n, position2), -Math.log(minorAlleleFractions.size()));
        }
    }
}