package org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter;


import org.apache.commons.math3.distribution.BinomialDistribution;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

public class ArtifactStatisticsScorer {
    
    final static double DEFAULT_BIASQP1=36;
    final static double DEFAULT_BIASQP2=1.5;

    /**
     *  Reduces the number of artifacts to cut based on the preAdapterQ score.
     *
     * @param preAdapterQ preAdapterQ score calculated by picard tool.  Must be larger than zero.  Phred score.
     * @param biasQP1 Inflection point for a sigmoid curve based on preAdapterQ.
     * @param biasQP2 If zero, then sharp cutoff at the preAdapterQ specified at biasQP1
     * @return multiplier to number of orientation bias artifact mutations to cut.  Always between 0 and 1.  Zero indicates that this sample probably has no artifacts, so nothing should be cut.
     */
    public static double calculateSuppressionFactorFromPreAdapterQ(final double preAdapterQ, final double biasQP1, final double biasQP2) {

        ParamUtils.isPositive(preAdapterQ, "preAdapter Q score must be positive and not zero.");
        ParamUtils.isPositiveOrZero(biasQP1, "bias Q shape parameters must be positive.");
        ParamUtils.isPositiveOrZero(biasQP2, "bias Q shape parameters must be positive.");

        // From matlab: fQ=1./(1+exp(biasQP2*(biasQ-biasQP1)));
        // biasQ is the preAdapterQ score
        // (biasQP2==0)&(biasQ>biasQP1) then return 0
        if (biasQP2 == 0) {
            return preAdapterQ > biasQP1 ? 0.0 : 1.0;
        } else {
            return 1.0/(1.0 + Math.exp(biasQP2*(preAdapterQ-biasQP1)));
        }
    }

    /**
     * See {@link #calculateSuppressionFactorFromPreAdapterQ(double, double, double) calculateSuppressionFactorFromPreAdapterQ}
     *
     *  Configured to behave exactly like the CGA Matlab OxoG filter in Firehose.  For those not internal to the Broad, this
     *   would be the configuration that is used in the vast majority of calls from the Broad, including TCGA.
     *
     * @param preAdapterQ See {@link #calculateSuppressionFactorFromPreAdapterQ(double, double, double) calculateSuppressionFactorFromPreAdapterQ}
     * @return See {@link #calculateSuppressionFactorFromPreAdapterQ(double, double, double) calculateSuppressionFactorFromPreAdapterQ}
     */
    public static double calculateSuppressionFactorFromPreAdapterQ(final double preAdapterQ) {
        return calculateSuppressionFactorFromPreAdapterQ(preAdapterQ, DEFAULT_BIASQP1, DEFAULT_BIASQP2);
    }

    /** p-value for being an artifact
     *
     * @param totalAltAlleleCount total number of alt reads
     * @param artifactAltAlleleCount alt read counts in a pair orientation that supports an artifact
     * @param biasP believed bias p value for a binomial distribution in artifact variants
     * @return p-value for this being an artifact
     */
    public static double calculateArtifactPValue(final int totalAltAlleleCount, final int artifactAltAlleleCount, final double biasP) {
        ParamUtils.isPositiveOrZero(biasP, "bias parameter must be positive or zero.");
        ParamUtils.isPositiveOrZero(totalAltAlleleCount, "total alt allele count must be positive or zero.");
        ParamUtils.isPositiveOrZero(artifactAltAlleleCount, "artifact supporting alt allele count must be positive or zero.");
        ParamUtils.isPositiveOrZero(totalAltAlleleCount-artifactAltAlleleCount, "Total alt count must be same or greater than the artifact alt count.");
        return new BinomialDistribution(null, totalAltAlleleCount, biasP).cumulativeProbability(artifactAltAlleleCount);
    }
}
