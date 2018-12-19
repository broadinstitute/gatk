package org.broadinstitute.hellbender.tools.walkers.mutect;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.utils.Dirichlet;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

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
     * @param log10Likelihoods matrix of alleles x reads
     * @param priorPseudocounts
     */
    public static double[] alleleFractionsPosterior(final RealMatrix log10Likelihoods, final double[] priorPseudocounts) {
        final int numberOfAlleles = log10Likelihoods.getRowDimension();
        Utils.validateArg(numberOfAlleles == priorPseudocounts.length, "Must have one pseudocount per allele.");

        double[] dirichletPosterior = new IndexRange(0, numberOfAlleles).mapToDouble(n -> 1.0);  // initialize flat posterior
        boolean converged = false;

        while(!converged) {
            // alleleCounts = \sum_r \bar{z}_r, where \bar{z}_r is an a-dimensional vector of the expectation of z_r with respect to q(f)
            final double[] alleleCounts = getEffectiveCounts(log10Likelihoods, dirichletPosterior);
            final double[] newDirichletPosterior = MathArrays.ebeAdd(alleleCounts, priorPseudocounts);
            converged = MathArrays.distance1(dirichletPosterior, newDirichletPosterior) < CONVERGENCE_THRESHOLD;
            dirichletPosterior = newDirichletPosterior;
        }

        return dirichletPosterior;
    }

    //same with flat prior
    public static double[] alleleFractionsPosterior(final RealMatrix log10Likelihoods) {
        final double[] flatPrior = new IndexRange(0, log10Likelihoods.getRowDimension()).mapToDouble(n -> 1);
        return alleleFractionsPosterior(log10Likelihoods, flatPrior);
    }


    /**
     * Given data log likelihoods and a Dirichlet prior for a categorical distribution, obtain the array of total
     * responsibilities for each category
     * @param log10Likelihoods
     * @param dirichletPrior
     * @return
     */
    @VisibleForTesting
    protected static double[] getEffectiveCounts(RealMatrix log10Likelihoods, double[] dirichletPrior) {
        final double[] effectiveLog10Weights = new Dirichlet(dirichletPrior).effectiveLog10MultinomialWeights();
        return MathUtils.sumArrayFunction(0, log10Likelihoods.getColumnDimension(),
                read -> MathUtils.posteriors(effectiveLog10Weights, log10Likelihoods.getColumn(read)));
    }

    /**
     * @param log10Likelihoods matrix of alleles x reads
     * @param priorPseudocounts
     */
    public static double log10Evidence(final RealMatrix log10Likelihoods, final double[] priorPseudocounts) {
        return log10Evidence(log10Likelihoods, priorPseudocounts, 0.0, -1);
    }


        /**
         * @param log10Likelihoods matrix of alleles x reads (NOTE: NON_REF allele is assumed to be last)
         * @param priorPseudocounts
         * @param alleleFractionThreshold lower bound of allele fractions to consider for non-ref likelihood
         */
    public static double log10Evidence(final RealMatrix log10Likelihoods, final double[] priorPseudocounts, final double alleleFractionThreshold, final int nonRefIndex) {
        final int numberOfAlleles = log10Likelihoods.getRowDimension();
        Utils.validateArg(numberOfAlleles == priorPseudocounts.length, "Must have one pseudocount per allele.");
        final double[] alleleFractionsPosterior = alleleFractionsPosterior(log10Likelihoods, priorPseudocounts);
        final double priorContribution = log10DirichletNormalization(priorPseudocounts);
        final double posteriorContribution = -log10DirichletNormalization(alleleFractionsPosterior);
        final double posteriorTotal = MathUtils.sum(alleleFractionsPosterior);
        double thresholdedPosteriorContribution = posteriorContribution;
        if (nonRefIndex > 0) {
            thresholdedPosteriorContribution += Math.log10(1-Beta.regularizedBeta(alleleFractionThreshold,
                    alleleFractionsPosterior[nonRefIndex], posteriorTotal - alleleFractionsPosterior[nonRefIndex]));
        }

        final double[] log10AlleleFractions = new Dirichlet(alleleFractionsPosterior).effectiveLog10MultinomialWeights();

        final double likelihoodsAndEntropyContribution = new IndexRange(0, log10Likelihoods.getColumnDimension()).sum(r -> {
            final double[] log10LikelihoodsForRead = log10Likelihoods.getColumn(r);
            final double[] responsibilities = MathUtils.posteriors(log10AlleleFractions, log10LikelihoodsForRead);
            final double entropyContribution = Arrays.stream(responsibilities).map(SomaticLikelihoodsEngine::xLog10x).sum();
            return likelihoodsContribution(log10LikelihoodsForRead, responsibilities) - entropyContribution;
        });

        return priorContribution + thresholdedPosteriorContribution + likelihoodsAndEntropyContribution;
    }

    private static double likelihoodsContribution(final double[] log10LikelihoodsForRead, final double[] responsibilities) {
        // this is a safe version of MathUtils.sum(MathArrays.ebeMultiply(log10LikelihoodsForRead, responsibilities))
        // in case the likelihood is zero, and the log likelihood in -Infinity, we have the responsibility is zero and the
        // contribution x * log(y), where x and y go to zero at the same rate (ie within a constant factor of each other
        // since the responsibility is related to the likelihood via the prior), is undefined but should be treated as zero.
        double result = 0;
        for (int n = 0; n < log10LikelihoodsForRead.length; n++) {
            result += (responsibilities[n] < NEGLIGIBLE_RESPONSIBILITY ? 0 : log10LikelihoodsForRead[n] * responsibilities[n]);
        }
        return result;
    }


    // same as above using the default flat prior
    public static double log10Evidence(final RealMatrix log10Likelihoods, final double minAF, final int nonRefIndex) {
        final double[] flatPrior = new IndexRange(0, log10Likelihoods.getRowDimension()).mapToDouble(n -> 1);
        return log10Evidence(log10Likelihoods, flatPrior, minAF, nonRefIndex);
    }

    private static double xLog10x(final double x) {
        return x < 1e-8 ? 0 : x * Math.log10(x);
    }

    public static double log10DirichletNormalization(final double... dirichletParams) {
        final double logNumerator = Gamma.logGamma(MathUtils.sum(dirichletParams));
        final double logDenominator = MathUtils.sum(MathUtils.applyToArray(dirichletParams, Gamma::logGamma));
        return MathUtils.logToLog10(logNumerator - logDenominator);
    }

}
