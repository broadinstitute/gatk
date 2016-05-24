package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.tools.exome.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Collection;
import java.util.function.Function;
import java.util.stream.IntStream;

import static java.lang.Math.log;
import static java.lang.Math.sqrt;
import static org.broadinstitute.hellbender.utils.MathUtils.log10Factorial;
import static org.broadinstitute.hellbender.utils.MathUtils.log10ToLog;

/**
 * Contains likelihood methods for the allele-fraction model.
 * See docs/CNVs/CNV-methods.pdf for a thorough description of the model.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AlleleFractionLikelihoods {
    private AlleleFractionLikelihoods() {}

    /**
     * Compute the log-likelihood of a alt reads and r ref reads given minor fraction f and gamma hyperparameters
     * (specifying the distribution on allelic biases) mu (mean) and beta (rate = mean/variance) and given
     * an alt minor, ref minor, or outlier indicator state.  Note that this is a partially collapsed log-likelihood in that the
     * latent variable corresponding to the allelic bias at this site has been marginalized out but the indicator
     * variable has not been marginalized out.
     * <p>
     * See docs/CNVs/CNV-methods.pdf for derivation.
     * <p>
     * Finally, note that this is a static method and does not get mu, beta, and minorFraction from an AlleleFractionState object
     * We need such functionality because MCMC evaluates the likelihood under proposed parameter changes.
     *
     * @param state allele fraction state
     * @param segment index of segment containijng this het site
     * @param count AllelicCount of alt and ref reads
     * @param indicator the hidden state (alt minor / ref minor / outlier)
     *
     * @return if indicator == ALT_MINOR:
     * <p>
     * log { [beta^alpha / Gamma(alpha)][(1-pi)/2] * int_{0 to infty} f^a * (1-f)^r * lambda^(alpha + r - 1) * exp(-beta*lambda)/(f + (1-f)*lambda)^n d lambda }
     * <p>
     * if indicator == REF_MINOR, same as ALT_MINOR but with f <--> 1 - f
     * <p>
     * if indicator == OUTLIER log {pi * a!r!/(n+1)!}
     * <p>
     * where alpha = mu*beta and n = a + r
     */
    public static double hetLogLikelihood(final AlleleFractionState state, final int segment, final AllelicCount count, final AlleleFractionIndicator indicator) {
        return hetLogLikelihood(state, segment, count, indicator, AllelicPanelOfNormals.EMPTY_PON);
    }

    /**
     * Computes {@link AlleleFractionLikelihoods#hetLogLikelihood} using allelic-bias priors derived from an
     * {@link AllelicPanelOfNormals}.  See docs/CNVs/CNV-methods.pdf for details.
     */
    public static double hetLogLikelihood(final AlleleFractionState state, final int segment, final AllelicCount count, final AlleleFractionIndicator indicator,
                                          final AllelicPanelOfNormals allelicPON) {
        final SimpleInterval site = count.getInterval();
        final double alpha;
        final double beta;
        if (allelicPON.equals(AllelicPanelOfNormals.EMPTY_PON)) {
            // if PON is not available, use alpha and beta learned from sample
            beta = state.meanBias() / state.biasVariance();
            alpha = state.meanBias() * beta;
        } else {
            // if PON is available, return alpha and beta at site if site is in PON, otherwise return MLE alpha and beta if site is not in PON
            alpha = allelicPON.getAlpha(site);
            beta = allelicPON.getBeta(site);
        }
        final double pi = state.outlierProbability();
        final double minorFraction = state.segmentMinorFraction(segment);
        final int a = count.getAltReadCount();
        final int r = count.getRefReadCount();

        return hetLogLikelihood(alpha, beta, minorFraction, pi, indicator, a, r);
    }

    private static double hetLogLikelihood(final double alpha, final double beta, final double minorFraction,
                                           final double outlierProbability, final AlleleFractionIndicator indicator,
                                           final int a, final int r) {
        if (indicator == AlleleFractionIndicator.OUTLIER) {
            return log(outlierProbability) + log10ToLog(log10Factorial(a) + log10Factorial(r) - log10Factorial(a + r + 1));
        } else {
            final double f = indicator == AlleleFractionIndicator.ALT_MINOR ? minorFraction : 1 - minorFraction;
            return log((1 - outlierProbability) / 2) + logIntegralOverAllelicBias(alpha, beta, f, a, r);
        }
    }

    /**
     * Calculates the log of the integral of the allelic-bias posterior (approximated as a Gamma distribution
     * with mode and curvature given by {@link AlleleFractionLikelihoods#biasPosteriorMode}
     * and {@link AlleleFractionLikelihoods#biasPosteriorCurvature}, respectively)
     * at given values of the hyperparameters for the allelic-bias Gamma-distribution prior,
     * the minor-allele fraction parameter, and the observed counts at a site.
     * See docs/CNVs/CNV-methods.pdf (where this quantity is referred to as phi) for details.
     * @param alpha alpha hyperparameter for allelic-bias Gamma-distribution prior
     * @param beta  beta hyperparameter for allelic-bias Gamma-distribution prior
     * @param f     minor-allele fraction
     * @param a     alt counts
     * @param r     ref counts
     */
    protected static double logIntegralOverAllelicBias(final double alpha, final double beta, final double f, final int a, final int r) {
        final double lambda0 = biasPosteriorMode(alpha, beta, f, a, r);
        final int n = a + r;
        final double kappa = biasPosteriorCurvature(alpha, f, r, n, lambda0);
        final double rho = biasPosteriorEffectiveAlpha(lambda0, kappa);
        final double tau = biasPosteriorEffectiveBeta(lambda0, kappa);
        final double logc = alpha*log(beta) - Gamma.logGamma(alpha) + a * log(f) + r * log(1 - f)
                + (r + alpha - rho) * log(lambda0) + (tau - beta) * lambda0 - n * log(f + (1 - f) * lambda0);
        return logc + Gamma.logGamma(rho) - rho * log(tau);
    }

    /**
     * Calculates the mode of the exact allelic-bias posterior at given values of the hyperparameters for the
     * * allelic-bias Gamma-distribution prior, the minor-allele fraction parameter, and the observed
     * counts at a site.  See docs/CNVs/CNV-methods.pdf (where this quantity is referred to as lambda_0) for details.
     * @param alpha alpha hyperparameter for allelic-bias Gamma-distribution prior
     * @param beta  beta hyperparameter for allelic-bias Gamma-distribution prior
     * @param f     minor-allele fraction
     * @param a     alt counts
     * @param r     ref counts
     */
    protected static double biasPosteriorMode(final double alpha, final double beta, final double f, final int a, final int r) {
        final double w = (1 - f) * (a - alpha + 1) + beta * f;
        return (sqrt(w * w + 4 * beta * f * (1 - f) * (r + alpha - 1)) - w) / (2 * beta * (1 - f));
    }

    /**
     * Calculates the curvature (second derivative at the mode) of the exact allelic-bias log posterior
     * at given values of the hyperparameters for the allelic-bias Gamma-distribution prior,
     * the minor-allele fraction parameter, and the observed counts at a site.
     * See docs/CNVs/CNV-methods.pdf (where this quantity is referred to as kappa) for details.
     * @param alpha     alpha hyperparameter for allelic-bias Gamma-distribution prior
     * @param f         minor-allele fraction
     * @param r         ref counts
     * @param n         total counts
     * @param lambda0   mode of allelic-bias posterior
     */
    protected static double biasPosteriorCurvature(final double alpha, final double f, final int r, final int n, final double lambda0) {
        final double y = (1 - f)/(f + (1 - f) * lambda0);
        return n * y * y - (r + alpha - 1) / (lambda0 * lambda0);
    }

    /**
     * Calculates the effective alpha hyperparameter for the Gamma-distribution approximation of the exact allelic-bias posterior.
     * See docs/CNVs/CNV-methods.pdf (where this quantity is referred to as rho) for details.
     * @param lambda0   mode of allelic-bias posterior
     * @param kappa     curvature of allelic-bias posterior
     */
    protected static double biasPosteriorEffectiveAlpha(final double lambda0, final double kappa) {
        return 1 - kappa * lambda0 * lambda0;
    }

    /**
     * Calculates the effective beta hyperparameter for the Gamma-distribution approximation of the exact allelic-bias posterior.
     * See docs/CNVs/CNV-methods.pdf (where this quantity is referred to as tau) for details.
     * @param lambda0   mode of allelic-bias posterior
     * @param kappa     curvature of allelic-bias posterior
     */
    protected static double biasPosteriorEffectiveBeta(final double lambda0, final double kappa) {
        return -kappa * lambda0;
    }

    /**
     * the log likelihood summed (marginalized) over indicator states, which we use in the fully collapsed model
     * in which latent variables (bias and indicator) are marginalized out
     *
     * @param state allele fraction state
     * @param segment index of segment containing this het site
     * @param count AllelicCount of alt and ref reads
     * @param allelicPON allelic panel of normals constructed from total alt and ref counts observed across all normals at each site
     * @return the log of the likelihood at this het site, marginalized over indicator states.
     */
    public static double collapsedHetLogLikelihood(final AlleleFractionState state, final int segment, final AllelicCount count, final AllelicPanelOfNormals allelicPON) {
        return GATKProtectedMathUtils.logSumExp(
                hetLogLikelihood(state, segment, count, AlleleFractionIndicator.ALT_MINOR, allelicPON),
                hetLogLikelihood(state, segment, count, AlleleFractionIndicator.REF_MINOR, allelicPON),
                hetLogLikelihood(state, segment, count, AlleleFractionIndicator.OUTLIER, allelicPON));
    }

    /**
     * the log-likelihood of all the hets in a segment
     *
     * @param state allele fraction state
     * @param segment index of segment containijng this het site
     * @param counts AllelicCount of alt and ref reads in this segment
     * @param allelicPON allelic panel of normals constructed from total alt and ref counts observed across all normals at each site
     * @return the sum of log-likelihoods over all het sites in a segment
     */
    public static double segmentLogLikelihood(final AlleleFractionState state, final int segment, final Collection<AllelicCount> counts, final AllelicPanelOfNormals allelicPON) {
        return counts.stream().mapToDouble(c -> collapsedHetLogLikelihood(state, segment, c, allelicPON)).sum();
    }

    /**
     * the total log likelihood of all segments
     * @param state current state
     * @param data data
     * @return sum of log likelihoods of all segments
     */
    public static double logLikelihood(final AlleleFractionState state, final AlleleFractionData data) {
        return IntStream.range(0, data.getNumSegments()).mapToDouble(s -> segmentLogLikelihood(state, s, data.getCountsInSegment(s), data.getPON())).sum();
    }

    protected static Function<Double, Double> segmentLogLikelihoodConditionalOnMinorFraction(final AlleleFractionState state,
                                                                                             final AlleleFractionData data, final int segment) {
        return minorFraction -> {
            final AlleleFractionState proposal = new AlleleFractionState(state.meanBias(), state.biasVariance(), state.outlierProbability(), minorFraction);
            return AlleleFractionLikelihoods.segmentLogLikelihood(proposal, 0, data.getCountsInSegment(segment), data.getPON());
        };
    }

    /**
     * Calculates the log-likelihood function for MLE initialization of the {@link AllelicPanelOfNormals}.
     * See docs/CNVs/CNV-methods.pdf for details.
     * @param meanBias      mean of the allelic-bias prior
     * @param biasVariance  variance of the allelic-bias prior
     * @param counts        counts at all sites in the {@link AllelicPanelOfNormals}
     */
    protected static double logLikelihoodForAllelicPanelOfNormals(final double meanBias, final double biasVariance, final AllelicCountCollection counts) {
        final double alpha = alpha(meanBias, biasVariance);
        final double beta = beta(meanBias, biasVariance);
        final double balancedMinorAlleleFraction = 0.5;
        return counts.getCounts().stream().mapToDouble(c -> AlleleFractionLikelihoods.logIntegralOverAllelicBias(alpha, beta, balancedMinorAlleleFraction, c.getAltReadCount(), c.getRefReadCount())).sum();
    }

    private static double alpha(final double meanBias, final double biasVariance) {
        return meanBias * meanBias / biasVariance;
    }

    private static double beta(final double meanBias, final double biasVariance) {
        return meanBias / biasVariance;
    }
}
