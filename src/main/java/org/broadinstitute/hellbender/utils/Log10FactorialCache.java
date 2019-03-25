package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.special.Gamma;

/**
 * Wrapper class so that the log10Factorial array is only calculated if it's used
 */
public final class Log10FactorialCache extends IntToDoubleFunctionCache {
    private static final int CACHE_SIZE = 10_000;

    @Override
    protected int maxSize() {
        return CACHE_SIZE;
    }

    @Override
    protected double compute(final int n) {
        return MathUtils.log10Gamma(n + 1);
    }
}
