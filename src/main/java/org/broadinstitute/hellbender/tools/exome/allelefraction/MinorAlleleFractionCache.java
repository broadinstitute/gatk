package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.broadinstitute.hellbender.tools.exome.AllelicCount;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * A helper class to maintain a cache of computed maximum-likelihood estimates of minor allele fraction
 * for {@link AllelicCount#estimateMinorAlleleFraction(double)}.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class MinorAlleleFractionCache {
    private MinorAlleleFractionCache() {}

    // constants for Brent optimization of minor allele fraction
    private static final MaxEval BRENT_MAX_EVAL = new MaxEval(100);
    private static final double MINOR_ALLELE_FRACTION_RELATIVE_TOLERANCE = 0.00001;
    private static final double MINOR_ALLELE_FRACTION_ABSOLUTE_TOLERANCE = 0.00001;
    private static final BrentOptimizer OPTIMIZER =
            new BrentOptimizer(MINOR_ALLELE_FRACTION_RELATIVE_TOLERANCE, MINOR_ALLELE_FRACTION_ABSOLUTE_TOLERANCE);

    private static final Map<MinorAlleleFractionCacheKey, Double> cache = new HashMap<>();

    public static double get(final int a, final int r, final double allelicBias) {
        final MinorAlleleFractionCacheKey key = new MinorAlleleFractionCacheKey(a, r, allelicBias);
        if (cache.containsKey(key)) {
            return cache.get(key);
        } else {
            final double minorAlleleFraction = estimateMinorAlleleFraction(a, r, allelicBias);
            cache.put(key, minorAlleleFraction);
            return minorAlleleFraction;
        }
    }

    //See docs/CNVs/CNV-methods.pdf for derivation of likelihood
    private static double estimateMinorAlleleFraction(final int a, final int r, final double bias) {
        final double altFraction = (double) a / (a + r);
        final double initialEstimate = altFraction < 0.5 ? altFraction : 1. - altFraction;
        final SearchInterval searchInterval = new SearchInterval(0.0, 0.5, initialEstimate);

        // work in log space to avoid underflow
        final UnivariateObjectiveFunction objective = new UnivariateObjectiveFunction(f -> {
            final double logf = Math.log(f);
            final double logOneMinusf = Math.log(1. - f);
            final double logBias = Math.log(bias);
            final double logAltMinorDenominatorTerm = Math.log(f + (1. - f)*bias);
            final double logRefMinorDenominatorTerm = Math.log(1. - f + f*bias);
            final double altMinorLogLikelihood =
                    a * logf + r * (logOneMinusf + logBias) - (a + r) * logAltMinorDenominatorTerm;
            final double refMinorLogLikelihood =
                    r * (logf + logBias) + a * logOneMinusf - (a + r) * logRefMinorDenominatorTerm;
            return GATKProtectedMathUtils.naturalLogSumExp(altMinorLogLikelihood, refMinorLogLikelihood);
        });

        return OPTIMIZER.optimize(objective, GoalType.MAXIMIZE, searchInterval, BRENT_MAX_EVAL).getPoint();
    }

    /**
     * A helper class to maintain a cache of computed maximum-likelihood estimates of minor allele fraction
     * for {@link AllelicCount#estimateMinorAlleleFraction(double)}.
     *
     * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
     * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
     */
    private static final class MinorAlleleFractionCacheKey {
        private final int a;
        private final int r;
        private final double bias;

        public MinorAlleleFractionCacheKey(final int a, final int r, final double bias) {
            this.a = a;
            this.r = r;
            this.bias = bias;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }

            final MinorAlleleFractionCacheKey that = (MinorAlleleFractionCacheKey) o;
            return a == that.a && r == that.r && Double.compare(that.bias, bias) == 0;
        }

        @Override
        public int hashCode() {
            final int[] hashes = new int[]{Integer.hashCode(a), Integer.hashCode(r), Double.hashCode(bias)};
            return Arrays.hashCode(hashes);
        }
    }
}