package org.broadinstitute.hellbender.tools.coveragemodel;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;

/**
 * Implements the {@link TargetLikelihoodCalculator} interface for the original XHMM-based germline model.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class XHMMTargetLikelihoodCalculator implements TargetLikelihoodCalculator {
    private final RandomGenerator rng;
    private final double emissionStandardDeviation;
    private final double deletionMean;
    private final double duplicationMean;
    private static final double NEUTRAL_MEAN = 0.0;

    public XHMMTargetLikelihoodCalculator(final double deletionMean, final double duplicationMean, final double emissionStdDev,
                                          final RandomGenerator rng) {
        this.rng = rng;
        this.emissionStandardDeviation = emissionStdDev;
        this.duplicationMean = duplicationMean;
        this.deletionMean = deletionMean;
    }

    @Override
    public double logLikelihood(final Target target, final double copyRatio, double coverageZScore) {
        return new NormalDistribution(rng, getEmissionMean(copyRatio), emissionStandardDeviation).logDensity(coverageZScore);
    }

    private double getEmissionMean(final double copyRatio) {
        if (copyRatio == CopyNumberTriState.NEUTRAL.copyRatio) {
            return NEUTRAL_MEAN;
        } else if (copyRatio == CopyNumberTriState.DUPLICATION.copyRatio) {
            return duplicationMean;
        } else if (copyRatio == CopyNumberTriState.DELETION.copyRatio){
            return deletionMean;
        } else {
            throw new IllegalArgumentException("The only valid copy ratios are those of CopyNumberTriState.");
        }
    }

    public double generateRandomZScoreData(final double copyRatio) {
        return getEmissionMean(copyRatio) + rng.nextGaussian() * emissionStandardDeviation;
    }

    public double deletionMean() { return deletionMean; }
    public double duplicationMean() { return duplicationMean; }
}
