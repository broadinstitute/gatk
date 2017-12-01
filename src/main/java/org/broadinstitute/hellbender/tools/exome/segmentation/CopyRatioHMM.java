package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.CauchyDistribution;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;
import java.util.stream.Collectors;

/**
 * The Copy Ratio Hidden Markov Model is a generative model describing latent CNV states
 * with different copy ratios.
 *
 * Hidden states are essentially copy ratios c, but for convenience we represent them
 * as integers 0, 1, . . . K - 1 and store copy ratios in a corresponding array.
 *
 * The model contains a memory length parameter representing the prior probability for the minor
 * allele fraction state to be "forgotten" as a function of distance d between consecutive SNPs --
 * the probability to remember a state is exp(-d/memoryLength).  If the state is forgotten, a new state
 * is chosen with probabilities given by an array of weights.
 *
 * Thus our transition probabilities are P(i -> j) = exp(-d/D) delta_{ij} + (1 - exp(-d/D)) weights[j]
 * where delta is the Kronecker delta and D is the memory length.
 *
 * Emission likelihoods as a function of copy ratio are assumed Gaussian in log2 space with a shared global variance.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class CopyRatioHMM extends ClusteringGenomicHMM<Double, Double> {
    private final double logCoverageCauchyWidth;
    private final List<AbstractRealDistribution> emissionDistributions;

    /**
     * @param log2CopyRatios array of log-2 copy ratios corresponding to the hidden states
     * @param memoryLength when consecutive SNPs are a distance d bases apart, the prior probability
     */
    public CopyRatioHMM(final List<Double> log2CopyRatios, final double memoryLength, final double logCoverageCauchyWidth) {
        super(log2CopyRatios, memoryLength);
        this.logCoverageCauchyWidth = logCoverageCauchyWidth;
        emissionDistributions = hiddenStates().stream()
                .map(n -> new CauchyDistribution(null, log2CopyRatios.get(n), logCoverageCauchyWidth)).collect(Collectors.toList());
    }

    /**
     * Visible for the segmenter
     */
    @Override
    public double logEmissionProbability(final Double data, final Integer state, final SimpleInterval position) {
        return emissionDistributions.get(state).logDensity(data);
    }

    /**
     * Visible for the segmenter
     */
    @Override
    public double logEmissionProbability(final Double data, final Double log2CopyRatio) {
        return logEmissionProbability(data, log2CopyRatio, logCoverageCauchyWidth);
    }

    /**
     * Visible for the segmenter
     */
    public static double logEmissionProbability(final Double data, final double log2CopyRatio, final double cauchyWidth) {
        return new CauchyDistribution(null, log2CopyRatio, cauchyWidth).logDensity(data);
    }
}