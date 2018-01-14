package org.broadinstitute.hellbender.tools.copynumber.models;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.IntStream;

import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.broadinstitute.hellbender.utils.MathUtils.log10Factorial;
import static org.broadinstitute.hellbender.utils.MathUtils.log10ToLog;

/**
 * Contains likelihood methods for the allele-fraction model.
 * See docs/CNVs/CNV-methods.pdf for a thorough description of the model.
 *
 * We can compute the log-likelihood of a alt reads and r ref reads given minor fraction f and gamma hyperparameters
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
 * <p>
 * if indicator == ALT_MINOR:
 * <p>
 * log { [beta^alpha / Gamma(alpha)][(1 - pi) / 2] * int_{0 to infty} f^a * (1 - f)^r * lambda^(alpha + r - 1) * exp(-beta * lambda)/(f + (1 - f) * lambda)^n d lambda }
 * <p>
 * if indicator == REF_MINOR, same as ALT_MINOR but with f <--> 1 - f
 * <p>
 * if indicator == OUTLIER log {pi * a!r!/(n+1)!}
 * <p>
 * where alpha = mu*beta and n = a + r.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class AlleleFractionLikelihoods {
    private static final double EPSILON = 1E-10;

    private static final FunctionCache logGammaCache = new FunctionCache(Gamma::logGamma);
    private static final FunctionCache logCache = new FunctionCache(FastMath::log);

    private static final class FunctionCache extends LinkedHashMap<Double, Double> {
        private static final long serialVersionUID = 19841647L;
        private static final int MAX_SIZE = 100_000;

        private final Function<Double, Double> mappingFunction;

        FunctionCache(final Function<Double, Double> mappingFunction) {
            this.mappingFunction = mappingFunction;
        }

        Double computeIfAbsent(final Double key) {
            return super.computeIfAbsent(key, mappingFunction);
        }

        @Override
        protected boolean removeEldestEntry(final Map.Entry<Double, Double> eldest) {
            return size() >= MAX_SIZE;
        }
    }

    private AlleleFractionLikelihoods() {}

    static double segmentLogLikelihood(final AlleleFractionGlobalParameters parameters,
                                       final double minorFraction,
                                       final List<AlleleFractionSegmentedData.IndexedAllelicCount> allelicCountsInSegment) {
        final double alpha = parameters.getAlpha();
        final double beta = parameters.getBeta();
        final double pi = parameters.getOutlierProbability();

        //we compute some quantities that will be reused
        final double logPi = logCache.computeIfAbsent(pi);
        final double logNotPi = logCache.computeIfAbsent((1 - pi) / 2);
        final double logcCommon = alpha * logCache.computeIfAbsent(beta) - logGammaCache.computeIfAbsent(alpha);
        final double majorFraction = 1 - minorFraction;
        final double logMinorFraction = log(minorFraction);
        final double logMajorFraction = log(majorFraction);

        double logLikelihood = 0.;
        for (final AllelicCount allelicCount : allelicCountsInSegment) {
            final int a = allelicCount.getAltReadCount();
            final int r = allelicCount.getRefReadCount();
            final int n = a + r;

            //alt-minor calculation
            final double lambda0AltMinor = biasPosteriorMode(alpha, beta, minorFraction, a, r);
            final double kappaAltMinor = biasPosteriorCurvature(alpha, minorFraction, r, n, lambda0AltMinor);
            final double rhoAltMinor = biasPosteriorEffectiveAlpha(lambda0AltMinor, kappaAltMinor);
            final double tauAltMinor = biasPosteriorEffectiveBeta(lambda0AltMinor, kappaAltMinor);
            final double logcAltMinor = logcCommon + a * logMinorFraction + r * logMajorFraction
                    + (r + alpha - rhoAltMinor) * log(lambda0AltMinor) + (tauAltMinor - beta) * lambda0AltMinor
                    - n * log(minorFraction + majorFraction * lambda0AltMinor);
            final double altMinorLogLikelihood = logNotPi + logcAltMinor + Gamma.logGamma(rhoAltMinor) - rhoAltMinor * log(tauAltMinor);

            //ref-minor calculation
            final double lambda0RefMinor = biasPosteriorMode(alpha, beta, majorFraction, a, r);
            final double kappaRefMinor = biasPosteriorCurvature(alpha, majorFraction, r, n, lambda0RefMinor);
            final double rhoRefMinor = biasPosteriorEffectiveAlpha(lambda0RefMinor, kappaRefMinor);
            final double tauRefMinor = biasPosteriorEffectiveBeta(lambda0RefMinor, kappaRefMinor);
            final double logcRefMinor = logcCommon + a * logMajorFraction + r * logMinorFraction
                    + (r + alpha - rhoRefMinor) * log(lambda0RefMinor) + (tauRefMinor - beta) * lambda0RefMinor
                    - n * log(majorFraction + minorFraction * lambda0RefMinor);
            final double refMinorLogLikelihood = logNotPi + logcRefMinor + Gamma.logGamma(rhoRefMinor) - rhoRefMinor * log(tauRefMinor);

            final double outlierLogLikelihood = logPi + log10ToLog(log10Factorial(a) + log10Factorial(r) - log10Factorial(a + r + 1));

            logLikelihood += MathUtils.logSumExp(altMinorLogLikelihood, refMinorLogLikelihood, outlierLogLikelihood);
        }
        return logLikelihood;
    }

    /**
     * The total log likelihood of all segments.
     */
    static double logLikelihood(final AlleleFractionGlobalParameters parameters,
                                final AlleleFractionState.MinorFractions minorFractions,
                                final AlleleFractionSegmentedData data) {
        return IntStream.range(0, data.getNumSegments())
                .mapToDouble(segment -> segmentLogLikelihood(parameters, minorFractions.get(segment), data.getIndexedAllelicCountsInSegment(segment)))
                .sum();
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
    private static double biasPosteriorMode(final double alpha, final double beta, final double f, final int a, final int r) {
        final double w = (1 - f) * (a - alpha + 1) + beta * f;
        return Math.max((sqrt(w * w + 4 * beta * f * (1 - f) * (r + alpha - 1)) - w) / (2 * beta * (1 - f)), EPSILON);
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
    private static double biasPosteriorCurvature(final double alpha, final double f, final int r, final int n, final double lambda0) {
        final double y = (1 - f) / (f + (1 - f) * lambda0);
        return n * y * y - (r + alpha - 1) / (lambda0 * lambda0);
    }

    /**
     * Calculates the effective alpha hyperparameter for the Gamma-distribution approximation of the exact allelic-bias posterior.
     * See docs/CNVs/CNV-methods.pdf (where this quantity is referred to as rho) for details.
     * @param lambda0   mode of allelic-bias posterior
     * @param kappa     curvature of allelic-bias posterior
     */
    private static double biasPosteriorEffectiveAlpha(final double lambda0, final double kappa) {
        return Math.max(1 - kappa * lambda0 * lambda0, EPSILON);
    }

    /**
     * Calculates the effective beta hyperparameter for the Gamma-distribution approximation of the exact allelic-bias posterior.
     * See docs/CNVs/CNV-methods.pdf (where this quantity is referred to as tau) for details.
     * @param lambda0   mode of allelic-bias posterior
     * @param kappa     curvature of allelic-bias posterior
     */
    private static double biasPosteriorEffectiveBeta(final double lambda0, final double kappa) {
        return Math.max(-kappa * lambda0, EPSILON);
    }

    private static double log(final double x) {
        return FastMath.log(Math.max(EPSILON, x));
    }
}