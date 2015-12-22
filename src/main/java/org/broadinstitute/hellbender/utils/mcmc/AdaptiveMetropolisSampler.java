package org.broadinstitute.hellbender.utils.mcmc;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.function.Function;

/**
 * Metropolis MCMC sampler using an adaptive step size that increases / decreases in order to decrease / increase acceptance
 * rate to some desired value.  (A general property of MCMC is that too-low acceptance rate are bad for obvious reasons
 * but too-high acceptance rates are also undesirable because it implies that steps are too small).
 * <p>
 * In order for the Markov chain to converge to the correct posterior distribution, adaptations to the step size must
 * vanish as sampling proceeds.
 * <p>
 * This sampling method is a very good black-box algorithm when we are reasonably confident that the sampled
 * conditional distribution is close to unimodal but otherwise unknown (i.e. not necessarily log-concave and with
 * unknown shape and width).
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AdaptiveMetropolisSampler {
    private static final double DEFAULT_OPTIMAL_ACCEPTANCE_RATE = 0.4;
    private static final double DEFAULT_TIME_SCALE = 20;
    private static final double DEFAULT_ADJUSTMENT_RATE = 1.0;
    private static final int RANDOM_SEED = 15;
    private static final RandomGenerator RNG = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

    private int iteration = 1;  // the amount of step size adjustment decreases with the iteration number
    private final double lowerBound;
    private final double upperBound;
    private double optimalAcceptanceRate = DEFAULT_OPTIMAL_ACCEPTANCE_RATE;

    //adjustments to the step size are scaled by adjustmentRate * timeScale / (timeScale + iteration)
    private final double adjustmentRate;
    private final double timeScale;

    private final AbstractRealDistribution normal = new NormalDistribution(RNG, 0, 1);
    private double stepSize;
    private double xCurrent;    //the current sampled value

    public AdaptiveMetropolisSampler(final double xInitial, final double initialStepSize, final double lowerBound,
            final double upperBound, final double adjustmentRate, final double timeScale) {
        xCurrent = xInitial;
        stepSize = initialStepSize;
        this.lowerBound = lowerBound;
        this.upperBound = upperBound;
        this.adjustmentRate = adjustmentRate;
        this.timeScale = timeScale;
    }

    public AdaptiveMetropolisSampler(final double xInitial, final double initialStepSize,
                                     final double lowerBound, final double upperBound) {
        this(xInitial, initialStepSize, lowerBound, upperBound, DEFAULT_ADJUSTMENT_RATE, DEFAULT_TIME_SCALE);
    }

    public double sample(final Function<Double, Double> logPDF) {
        final double proposal = xCurrent + stepSize * normal.sample();
        final double acceptanceProbability = (proposal < lowerBound || upperBound < proposal) ? 0
                : Math.min(1, Math.exp(logPDF.apply(proposal) - logPDF.apply(xCurrent)));

        //adjust stepSize larger/smaller to decrease/increase the acceptance rate
        final double correctionFactor = (acceptanceProbability - optimalAcceptanceRate) * adjustmentRate * (timeScale / (timeScale + iteration));
        stepSize *= Math.exp(correctionFactor);
        iteration++;
        return RNG.nextDouble() < acceptanceProbability ? proposal : xCurrent;
    }

    /**
     * Generate multiple samples from the probability density function.
     * @param numSamples    number of samples to generate
     * @param numBurnIn    number of samples to discard
     * @return              samples drawn from the probability density function
     */
    public List<Double> sample(final Function<Double, Double> logPDF, final int numSamples, final int numBurnIn) {
        final List<Double> samples = new ArrayList<>(numSamples);
        for (int i = 0; i < numSamples; i++) {
            xCurrent = sample(logPDF);
            if (i > numBurnIn) {
                samples.add(xCurrent);
            }
        }
        return samples;
    }
}