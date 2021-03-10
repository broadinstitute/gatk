package org.broadinstitute.hellbender.tools.walkers.mutect;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.utils.*;

import java.util.Arrays;

/**
 * Created by David Benjamin on 3/9/17.
 */
public class SomaticLikelihoodsEngine {

    public static final double CONVERGENCE_THRESHOLD = 0.001;
    private static double NEGLIGIBLE_RESPONSIBILITY = 1.0e-10;

    /**
     * Given a likelihoods matrix, calculate the parameters of the Dirichlet posterior distribution on their allele
     * fractions, which define a discrete distribution.
     * @param logLikelihoods matrix of alleles x reads
     * @param priorPseudocounts
     */
    public static double[] alleleFractionsPosterior(final RealMatrix logLikelihoods, final double[] priorPseudocounts, final double[] weights) {
        final int numberOfAlleles = logLikelihoods.getRowDimension();
        Utils.validateArg(numberOfAlleles == priorPseudocounts.length, "Must have one pseudocount per allele.");

        double[] dirichletPosterior = new IndexRange(0, numberOfAlleles).mapToDouble(n -> 1.0);  // initialize flat posterior
        boolean converged = false;

        while(!converged) {
            // alleleCounts = \sum_r \bar{z}_r, where \bar{z}_r is an a-dimensional vector of the expectation of z_r with respect to q(f)
            final double[] alleleCounts = getEffectiveCounts(logLikelihoods, dirichletPosterior, weights);
            final double[] newDirichletPosterior = MathArrays.ebeAdd(alleleCounts, priorPseudocounts);
            converged = MathArrays.distance1(dirichletPosterior, newDirichletPosterior)/MathUtils.sum(newDirichletPosterior) < CONVERGENCE_THRESHOLD;
            dirichletPosterior = newDirichletPosterior;
        }

        return dirichletPosterior;
    }

    public static double[] alleleFractionsPosterior(final RealMatrix logLikelihoods, final double[] priorPseudocounts) {
        return alleleFractionsPosterior(logLikelihoods, priorPseudocounts, null);
    }


        /**
         * Given data log likelihoods and a Dirichlet prior for a categorical distribution, obtain the array of total
         * responsibilities for each category
         * @param logLikelihoods
         * @param dirichletPrior
         * @return
         */
    @VisibleForTesting
    protected static double[] getEffectiveCounts(RealMatrix logLikelihoods, double[] dirichletPrior, final double[] weights) {
        final double[] effectiveLogWeights = new Dirichlet(dirichletPrior).effectiveLogMultinomialWeights();
        return MathUtils.sumArrayFunction(0, logLikelihoods.getColumnDimension(),
                read -> {
                    final double[] unweighted = NaturalLogUtils.posteriors(effectiveLogWeights, logLikelihoods.getColumn(read));
                    return weights == null ? unweighted : MathUtils.applyToArrayInPlace(unweighted, d -> d * weights[read]);
                });
    }

    protected static double[] getEffectiveCounts(RealMatrix logLikelihoods, double[] dirichletPrior) {
        return getEffectiveCounts(logLikelihoods, dirichletPrior, null);
    }


        /**
         * @param logLikelihoods matrix of alleles x reads (NOTE: NON_REF allele is assumed to be last)
         * @param priorPseudocounts
         */
    public static double logEvidence(final RealMatrix logLikelihoods, final double[] priorPseudocounts) {
        final int numberOfAlleles = logLikelihoods.getRowDimension();
        Utils.validateArg(numberOfAlleles == priorPseudocounts.length, "Must have one pseudocount per allele.");
        final double[] alleleFractionsPosterior = alleleFractionsPosterior(logLikelihoods, priorPseudocounts);
        final double priorContribution = logDirichletNormalization(priorPseudocounts);
        final double posteriorContribution = -logDirichletNormalization(alleleFractionsPosterior);

        final double[] logAlleleFractions = new Dirichlet(alleleFractionsPosterior).effectiveLogMultinomialWeights();

        final double likelihoodsAndEntropyContribution = new IndexRange(0, logLikelihoods.getColumnDimension()).sum(r -> {
            final double[] logLikelihoodsForRead = logLikelihoods.getColumn(r);
            final double[] responsibilities = NaturalLogUtils.posteriors(logAlleleFractions, logLikelihoodsForRead);
            final double entropyContribution = Arrays.stream(responsibilities).map(SomaticLikelihoodsEngine::xLogx).sum();
            return likelihoodsContribution(logLikelihoodsForRead, responsibilities) - entropyContribution;
        });

        return priorContribution + posteriorContribution + likelihoodsAndEntropyContribution;
    }

    private static double likelihoodsContribution(final double[] logLikelihoodsForRead, final double[] responsibilities) {
        // this is a safe version of MathUtils.sum(MathArrays.ebeMultiply(logLikelihoodsForRead, responsibilities))
        // in case the likelihood is zero, and the log likelihood in -Infinity, we have the responsibility is zero and the
        // contribution x * log(y), where x and y go to zero at the same rate (ie within a constant factor of each other
        // since the responsibility is related to the likelihood via the prior), is undefined but should be treated as zero.
        double result = 0;
        for (int n = 0; n < logLikelihoodsForRead.length; n++) {
            result += (responsibilities[n] < NEGLIGIBLE_RESPONSIBILITY ? 0 : logLikelihoodsForRead[n] * responsibilities[n]);
        }
        return result;
    }

    private static double xLogx(final double x) {
        return x < 1e-8 ? 0 : x * Math.log(x);
    }

    public static double logDirichletNormalization(final double... dirichletParams) {
        final double logNumerator = Gamma.logGamma(MathUtils.sum(dirichletParams));
        final double logDenominator = MathUtils.sum(MathUtils.applyToArray(dirichletParams, Gamma::logGamma));
        return logNumerator - logDenominator;
    }

}
