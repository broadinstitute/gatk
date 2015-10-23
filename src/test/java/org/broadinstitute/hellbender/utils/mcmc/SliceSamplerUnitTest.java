package org.broadinstitute.hellbender.utils.mcmc;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.List;
import java.util.Random;
import java.util.function.Function;


/**
 * Unit tests for {@link SliceSampler}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SliceSamplerUnitTest {
    private static final int RANDOM_SEED = 42;
    private static final RandomGenerator rng =
            RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

    private static double relativeError(final double x, final double xTrue) {
        return Math.abs((x - xTrue) / xTrue);
    }

    /**
     * Test slice sampling of a normal distribution.  Checks that input mean and standard deviation are recovered
     * by 20000 samples to a relative error of 0.05% and 2%, respectively.
     */
    @Test
    public void testSliceSamplingOfNormalDistribution() {
        final double mean = 5.;
        final double standardDeviation = 0.75;
        final NormalDistribution normalDistribution = new NormalDistribution(mean, standardDeviation);
        final Function<Double, Double> normalLogPDF = normalDistribution::logDensity;

        final double xInitial = 1.;
        final double xMin = Double.NEGATIVE_INFINITY;
        final double xMax = Double.POSITIVE_INFINITY;
        final double width = 0.5;
        final int numSamples = 20000;
        final SliceSampler normalSampler = new SliceSampler(rng, normalLogPDF, xInitial, xMin, xMax, width);
        final List<Double> samples = normalSampler.sample(numSamples);

        final double sampleMean = new Mean().evaluate(Doubles.toArray(samples));
        final double sampleStandardDeviation = new StandardDeviation().evaluate(Doubles.toArray(samples));
        Assert.assertEquals(relativeError(sampleMean, mean), 0., 0.0005);
        Assert.assertEquals(relativeError(sampleStandardDeviation, standardDeviation), 0., 0.02);
    }
}