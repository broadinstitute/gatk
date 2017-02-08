package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionGlobalParameters;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionLikelihoods;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.List;

/**
 * The Allele Fraction Hidden Markov Model is a generative model describing latent CNV states
 * with different minor allele fractions and observed ref and alt allele counts at a set of
 * het sites.
 *
 * Hidden states are essentially minor allele fractions f, but for convenience we represent them
 * as integers 0, 1, . . . K - 1 and store minor allele fractions in a corresponding array.
 *
 * The model contains a memory length parameter representing the prior probability for the minor
 * allele fraction state to be "forgotten" as a function of distance d between consecutive SNPs --
 * the probability to remember a state is exp(-d/memoryLength).  If the state is forgotten, a new state
 * is chosen with probabilities given by an array of weights.
 *
 * Thus our transition probabilities are P(i -> j) = exp(-d/D) delta_{ij} + (1 - exp(-d/D)) / K
 * where delta is the Kronecker delta and D is the memory length.
 *
 * Emission likelihoods as a function of minor allele fraction are defined and described in
 * {@link AlleleFractionLikelihoods}.  These likelihoods require this class to contain a set of
 * {@link AlleleFractionGlobalParameters} and an {@link AllelicPanelOfNormals}.
 *
 * This model, along with the likelihoods model, is described thoroughly in docs/CNVs/CNV-methods.pdf
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionHMM extends ClusteringGenomicHMM<AllelicCount, Double> {
    public static final double DEFAULT_OUTLIER_PROBABILITY = 0.01; //TODO: make this a variable, perhaps, but maybe not because it's just a garbage absorber
    private static final double EPSILON = 1E-10;
    private final double log10OutlierProbability;
    private final double log10NonOutlierProbability;

    /**
     * @param minorAlleleFractions array of minor allele fractions corresponding to the hidden states
     * @param memoryLength when consecutive SNPs are a distance d bases apart, the prior probability
     *                     for memory of the CNV state to be kept is exp(-d/memoryLength)
     */
    public AlleleFractionHMM(final List<Double> minorAlleleFractions, final double memoryLength, final double outlierProbability) {
        super(minorAlleleFractions, memoryLength);
        minorAlleleFractions.forEach(f -> ParamUtils.inRange(f, 0, 0.5, "minor fractions must be between 0 and 1/2, found " + f));
        log10OutlierProbability = Math.log10(outlierProbability);
        log10NonOutlierProbability = Math.log10(1 - outlierProbability);
    }

    public AlleleFractionHMM(final List<Double> minorAlleleFractions, final double memoryLength) {
        this(minorAlleleFractions, memoryLength, DEFAULT_OUTLIER_PROBABILITY);
    }

    @Override
    public double logEmissionProbability(final AllelicCount data, final Integer state, final SimpleInterval position) {
        return logEmissionProbability(data, getHiddenStateValue(state));
    }

    /**
     * Visible for {@link AlleleFractionSegmenter}
     */
    @Override
    public double logEmissionProbability(final AllelicCount data, final Double minorFraction) {
        return logEmissionProbability(data, minorFraction, log10OutlierProbability, log10NonOutlierProbability);
    }

    public static double logEmissionProbability(final AllelicCount data, final Double minorFraction,
                                                final double log10OutlierProbability, final double log10NonOutlierProbability) {
        final int altCount = data.getAltReadCount();
        final int refCount = data.getRefReadCount();
        final int totalCount = altCount + refCount;
        final double log10MinorFraction = Math.log10(Math.max(minorFraction, EPSILON));
        final double log10AltMinorLikelihood = MathUtils.log10BinomialProbability(totalCount, altCount, log10MinorFraction)
                + log10NonOutlierProbability + MathUtils.LOG10_ONE_HALF;
        final double log10RefMinorLikelihood = MathUtils.log10BinomialProbability(totalCount, refCount, log10MinorFraction)
                + log10NonOutlierProbability + MathUtils.LOG10_ONE_HALF;
        final double log10OutlierLikelihood = log10OutlierProbability - MathUtils.LOG10_ONE_HALF;   //flat outlier distribution over [0, 0.5]

        return MathUtils.log10ToLog(MathUtils.log10SumLog10(new double[] {log10AltMinorLikelihood, log10RefMinorLikelihood, log10OutlierLikelihood}));
    }

    public double getMinorAlleleFraction(final int state) { return getHiddenStateValue(state); }

}
