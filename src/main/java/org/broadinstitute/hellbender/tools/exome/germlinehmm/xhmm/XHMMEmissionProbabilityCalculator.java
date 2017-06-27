package org.broadinstitute.hellbender.tools.exome.germlinehmm.xhmm;

import com.beust.jcommander.internal.Nullable;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.broadinstitute.barclay.utils.Utils;
import org.broadinstitute.hellbender.tools.coveragemodel.interfaces.TargetLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.germlinehmm.CopyNumberTriState;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import java.util.Random;

/**
 * Implements the {@link TargetLikelihoodCalculator} interface for the original XHMM-based germline model.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class XHMMEmissionProbabilityCalculator implements TargetLikelihoodCalculator<XHMMEmissionData> {
    private final double emissionStandardDeviation;
    private final double deletionMean;
    private final double duplicationMean;
    private static final double NEUTRAL_MEAN = 0.0;

    public XHMMEmissionProbabilityCalculator(final double deletionMean, final double duplicationMean, final double emissionStdDev) {
        this.deletionMean = ParamUtils.isNegativeOrZero(deletionMean, "Deletion coverage shift must be negative.");
        this.duplicationMean = ParamUtils.isPositiveOrZero(duplicationMean, "Duplication coverage shift must be positive");
        emissionStandardDeviation = ParamUtils.isPositive(emissionStdDev, "Emission standard deviation must be positive");
    }

    @Override
    public double logLikelihood(@Nonnull final XHMMEmissionData emissionData, final double copyRatio, final Target target) {
        return new NormalDistribution(null, getEmissionMean(copyRatio), emissionStandardDeviation)
                .logDensity(emissionData.getCoverageZScore());
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

    public double generateRandomZScoreData(final double copyRatio, final Random rng) {
        Utils.nonNull(rng, "Random z-score sampling is called by the random number generator is null. This" +
                " instance of the class can not produce random samples.");
        return getEmissionMean(copyRatio) + rng.nextGaussian() * emissionStandardDeviation;
    }

    public double deletionMean() { return deletionMean; }
    public double duplicationMean() { return duplicationMean; }
}
