package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.special.Beta;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.List;
import java.util.Optional;

public class PolymeraseSlippageFilter extends Mutect2VariantFilter {
    private final int minSlippageLength;
    private final double slippageRate;


    public PolymeraseSlippageFilter(final int minSlippageLength, final double slippageRate) {
        this.minSlippageLength = minSlippageLength;
        this.slippageRate = slippageRate;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public double calculateErrorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringEngine) {

        final int[] rpa = vc.getAttributeAsList(GATKVCFConstants.REPEATS_PER_ALLELE_KEY).stream()
                .mapToInt(o -> Integer.parseInt(String.valueOf(o))).toArray();
        if (rpa.length < 2) {
            return 0;
        }
        final String ru = vc.getAttributeAsString(GATKVCFConstants.REPEAT_UNIT_KEY, "");

        final int referenceSTRBaseCount = ru.length() * rpa[0];
        final int numPCRSlips = rpa[0] - rpa[1];
        if (referenceSTRBaseCount >= minSlippageLength && Math.abs(numPCRSlips) == 1) {
            // calculate the p-value that out of n reads we would have at least k slippage reads
            // if this p-value is small we keep the variant (reject the PCR slippage hypothesis)
            final int[] ADs = filteringEngine.sumADsOverSamples(vc, true, false);
            if (ADs == null || ADs.length < 2) {
                return 0;
            }
            final int depth = (int) MathUtils.sum(ADs);

            final int altCount = (int) MathUtils.sum(ADs) - ADs[0];
            final double log10SomaticLikelihood = filteringEngine.getSomaticClusteringModel().log10LikelihoodGivenSomatic(depth, altCount);

            double likelihoodGivenSlippageArtifact;
            try {
                likelihoodGivenSlippageArtifact = Beta.regularizedBeta(slippageRate, ADs[1] + 1, ADs[0] + 1);
            } catch (final MaxCountExceededException e) {
                //if the special function can't be computed, use a binomial with fixed probability
                likelihoodGivenSlippageArtifact = new BinomialDistribution(null, depth, slippageRate).probability(ADs[1]);
            }

            final double log10Odds = log10SomaticLikelihood - Math.log10(likelihoodGivenSlippageArtifact);
            return filteringEngine.posteriorProbabilityOfError(vc, log10Odds, 0);
        } else {
            return 0;
        }
    }

    @Override
    public Optional<String> phredScaledPosteriorAnnotationName() {
        return Optional.of(GATKVCFConstants.POLYMERASE_SLIPPAGE_QUAL_VCF_ATTRIBUTE);
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.POLYMERASE_SLIPPAGE;
    }

    @Override
    protected List<String> requiredAnnotations() {
        return Arrays.asList(GATKVCFConstants.REPEATS_PER_ALLELE_KEY, GATKVCFConstants.REPEAT_UNIT_KEY);
    }
}
