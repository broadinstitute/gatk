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
            final double likelihoodsContribution = MathUtils.sum(MathArrays.ebeMultiply(log10LikelihoodsForRead, responsibilities));
            final double entropyContribution = Arrays.stream(responsibilities).map(SomaticLikelihoodsEngine::xLog10x).sum();
            return likelihoodsContribution - entropyContribution;
        });

        return priorContribution + thresholdedPosteriorContribution + likelihoodsAndEntropyContribution;
    }

    // same as above using the default flat prior
    public static double log10Evidence(final RealMatrix log10Likelihoods, final double minAF, final int nonRefIndex) {
        final double[] flatPrior = new IndexRange(0, log10Likelihoods.getRowDimension()).mapToDouble(n -> 1);
        return log10Evidence(log10Likelihoods, flatPrior, minAF, nonRefIndex);
    }

    private static double xLog10x(final double x) {
        return x < 1e-8 ? 0 : x * Math.log10(x);
    }

    public static double log10DirichletNormalization(final double[] dirichletParams) {
        final double logNumerator = Gamma.logGamma(MathUtils.sum(dirichletParams));
        final double logDenominator = MathUtils.sum(MathUtils.applyToArray(dirichletParams, Gamma::logGamma));
        return MathUtils.logToLog10(logNumerator - logDenominator);
    }

}
