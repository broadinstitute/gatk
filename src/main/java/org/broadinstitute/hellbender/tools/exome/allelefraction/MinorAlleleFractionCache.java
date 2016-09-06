package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.OptimizationUtils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.function.Function;

/**
 * A helper class to maintain a cache of computed maximum-likelihood estimates of minor allele fraction
 * for {@link AllelicCount#estimateMinorAlleleFraction(double)}.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class MinorAlleleFractionCache {
    private MinorAlleleFractionCache() {}

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

        // work in log space to avoid underflow
        final Function<Double, Double> objective = f -> {
            final double logf = Math.log(f);
            final double logOneMinusf = Math.log(1. - f);
            final double logBias = Math.log(bias);
            final double logAltMinorDenominatorTerm = Math.log(f + (1. - f)*bias);
            final double logRefMinorDenominatorTerm = Math.log(1. - f + f*bias);
            final double altMinorLogLikelihood =
                    a * logf + r * (logOneMinusf + logBias) - (a + r) * logAltMinorDenominatorTerm;
            final double refMinorLogLikelihood =
                    r * (logf + logBias) + a * logOneMinusf - (a + r) * logRefMinorDenominatorTerm;
            return GATKProtectedMathUtils.logSumExp(altMinorLogLikelihood, refMinorLogLikelihood);
        };
        try {
            return OptimizationUtils.argmax(objective, AlleleFractionState.MIN_MINOR_FRACTION, AlleleFractionState.MAX_MINOR_FRACTION, initialEstimate);
        } catch (final TooManyEvaluationsException ex){
            return initialEstimate;
        }
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