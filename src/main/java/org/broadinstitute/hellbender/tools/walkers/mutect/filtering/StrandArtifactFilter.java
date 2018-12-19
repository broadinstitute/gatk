package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.BetaBinomialDistribution;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.Collectors;

public class StrandArtifactFilter extends Mutect2VariantFilter {
    // beta prior on strand bias allele fraction
    private double INITIAL_ALPHA_STRAND = 1.0;
    private double INITIAL_BETA_STRAND = 6.0;

    private double alphaStrand = INITIAL_ALPHA_STRAND;
    private double betaStrand = INITIAL_BETA_STRAND;

    // beta prior on sequencing error allele fraction
    private static final double ALPHA_SEQ = 10;
    private static final double BETA_SEQ = 2000;

    private static final double INITIAL_STRAND_ARTIFACT_PRIOR = 0.01;

    private double strandArtifactPrior = INITIAL_STRAND_ARTIFACT_PRIOR;

    private static final double ARTIFACT_PSEUDOCOUNT = 1;

    private static final double NON_ARTIFACT_PSEUDOCOUNT = 10;

    private final List<EStep> eSteps = new ArrayList<>();

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public double calculateErrorProbability(final VariantContext vc, final Mutect2FilteringEngine filteringEngine) {
        final EStep probabilities = calculateArtifactProbabilities(vc, filteringEngine);
        return probabilities.forwardArtifactProbability + probabilities.reverseArtifactProbability;
    }

    public EStep calculateArtifactProbabilities(final VariantContext vc, final Mutect2FilteringEngine filteringEngine) {
        final int[] totalCounts = filteringEngine.sumADsOverSamples(vc, true, true);
        final int[] forwardCounts = vc.getAttributeAsIntList(GATKVCFConstants.FORWARD_STRAND_COUNT_KEY, 0).stream().mapToInt(x->x).toArray();
        final int forwardCount = (int) MathUtils.sum(forwardCounts);
        final int reverseCount = (int) MathUtils.sum(totalCounts) - forwardCount;
        final int forwardAltCount = forwardCount - forwardCounts[0];
        final int reverseAltCount = (int) MathUtils.sum(totalCounts) - totalCounts[0] - forwardAltCount;

        return strandArtifactProbability(strandArtifactPrior, forwardCount, reverseCount, forwardAltCount, reverseAltCount);

    }

    @Override
    protected void accumulateDataForLearning(final VariantContext vc, final ErrorProbabilities errorProbabilities, final Mutect2FilteringEngine filteringEngine) {
        final EStep eStep = calculateArtifactProbabilities(vc, filteringEngine);
        eSteps.add(eStep);
    }

    @Override
    protected void clearAccumulatedData() {
        eSteps.clear();
    }

    @Override
    protected void learnParameters() {
        final List<EStep> potentialArtifacts = eSteps.stream()
                .filter(eStep -> eStep.getArtifactProbability() > 0.1).collect(Collectors.toList());
        final double totalArtifacts = potentialArtifacts.stream().mapToDouble(EStep::getArtifactProbability).sum();
        final double totalNonArtifacts = eSteps.stream().mapToDouble(e -> 1 - e.getArtifactProbability()).sum();
        strandArtifactPrior = (totalArtifacts + ARTIFACT_PSEUDOCOUNT) / (totalNonArtifacts + NON_ARTIFACT_PSEUDOCOUNT);


        final double artifactAltCount = potentialArtifacts.stream()
                .mapToDouble(e -> e.forwardArtifactProbability * e.forwardAltCount + e.reverseArtifactProbability * e.reverseAltCount)
                .sum();

        final double artifactDepth = potentialArtifacts.stream()
                .mapToDouble(e -> e.forwardArtifactProbability * e.forwardCount + e.reverseArtifactProbability * e.reverseCount)
                .sum();

        final double artifactBetaMean = (artifactAltCount + INITIAL_ALPHA_STRAND) / (artifactDepth + INITIAL_ALPHA_STRAND + INITIAL_BETA_STRAND);

        // We do the M step for the beta prior on the artifact allele fraction by brute force single-parameter optimization.
        // By estimating the mean empirically as above we can fix mean = alpha / (alpha + beta), hence beta = (1/mean - 1) * alpha.
        // This lets us do single-parameter optimization on alpha with beta/alpha fixed.
        // brute force optimization is fairly cheap because the objective includes only calls that show some evidence of strand bias.
        final DoubleUnaryOperator objective = alpha -> {
            final double beta = (1 / artifactBetaMean - 1) * alpha;
            return potentialArtifacts.stream()
                    .mapToDouble( e -> e.getForwardArtifactProbability() * artifactStrandLogLikelihood(e.forwardCount, e.forwardAltCount, alpha, beta)
                            + e.getReverseArtifactProbability() * artifactStrandLogLikelihood(e.reverseCount, e.reverseAltCount, alpha, beta))
                    .sum();
        };

        alphaStrand = OptimizationUtils.max(objective, 0.01, 100, INITIAL_ALPHA_STRAND, 0.01, 0.01, 100).getPoint();
        betaStrand = (1/artifactBetaMean - 1)*alphaStrand;
        // free up memory
        eSteps.clear();
    }

    @VisibleForTesting
    EStep strandArtifactProbability(final double strandArtifactPrior, int forwardCount, int reverseCount, int forwardAltCount, int reverseAltCount) {
        final double forwardLogLikelihood = artifactStrandLogLikelihood(forwardCount, forwardAltCount)
                + nonArtifactStrandLogLikelihood(reverseCount, reverseAltCount);
        final double reverseLogLikelihood = artifactStrandLogLikelihood(reverseCount, reverseAltCount)
                + nonArtifactStrandLogLikelihood(forwardCount, forwardAltCount);
        final double noneLogLikelihood = CombinatoricsUtils.binomialCoefficientLog(forwardCount, forwardAltCount)
                + CombinatoricsUtils.binomialCoefficientLog(reverseCount, reverseAltCount)
                - CombinatoricsUtils.binomialCoefficientLog(forwardCount + reverseCount, forwardAltCount + reverseAltCount)
                + new BetaBinomialDistribution(null, 1, 1, forwardCount + reverseCount).logProbability(forwardAltCount + reverseAltCount);

        final double forwardLogPrior = Math.log(strandArtifactPrior/2);
        final double reverseLogPrior = Math.log(strandArtifactPrior/2);
        final double noneLogPrior = Math.log(1 - strandArtifactPrior);

        final double[] forwardReverseNoneProbs = MathUtils.normalizeLog10(new double[] {(forwardLogLikelihood + forwardLogPrior) * MathUtils.LOG10_OF_E,
                (reverseLogLikelihood + reverseLogPrior)* MathUtils.LOG10_OF_E, (noneLogLikelihood + noneLogPrior) * MathUtils.LOG10_OF_E}, false, true);

        return new EStep(forwardReverseNoneProbs[0], forwardReverseNoneProbs[1], forwardCount, reverseCount, forwardAltCount, reverseAltCount);
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() {
        return Arrays.asList(GATKVCFConstants.FORWARD_STRAND_COUNT_KEY);
    }

    private double artifactStrandLogLikelihood(final int strandCount, final int strandAltCount) {
        return artifactStrandLogLikelihood(strandCount, strandAltCount, alphaStrand, betaStrand);
    }

    private static double artifactStrandLogLikelihood(final int strandCount, final int strandAltCount, final double alpha, final double beta) {
        return new BetaBinomialDistribution(null, alpha, beta, strandCount).logProbability(strandAltCount);
    }

    private double nonArtifactStrandLogLikelihood(final int strandCount, final int strandAltCount) {
        return new BetaBinomialDistribution(null, ALPHA_SEQ, BETA_SEQ, strandCount).logProbability(strandAltCount);
    }

    @Override
    public Optional<String> phredScaledPosteriorAnnotationName() {
        return Optional.of(GATKVCFConstants.STRAND_QUAL_VCF_ATTRIBUTE);
    }

    public static final class EStep {
        private double forwardArtifactProbability;
        private double reverseArtifactProbability;
        private int forwardCount;
        private int reverseCount;
        private int forwardAltCount;
        private int reverseAltCount;

        public EStep(double forwardArtifactProbability, double reverseArtifactProbability, int forwardCount, int reverseCount, int forwardAltCount, int reverseAltCount) {
            this.forwardArtifactProbability = forwardArtifactProbability;
            this.reverseArtifactProbability = reverseArtifactProbability;
            this.forwardCount = forwardCount;
            this.reverseCount = reverseCount;
            this.forwardAltCount = forwardAltCount;
            this.reverseAltCount = reverseAltCount;
        }

        public double getForwardArtifactProbability() {
            return forwardArtifactProbability;
        }

        public double getReverseArtifactProbability() {
            return reverseArtifactProbability;
        }

        public double getArtifactProbability() {
            return getForwardArtifactProbability() + getReverseArtifactProbability();
        }

        public int getForwardCount() {
            return forwardCount;
        }

        public int getReverseCount() {
            return reverseCount;
        }

        public int getForwardAltCount() {
            return forwardAltCount;
        }

        public int getReverseAltCount() {
            return reverseAltCount;
        }
    }

}
