package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.special.Gamma;

public final class DigammaCache extends IntToDoubleFunctionCache {
    private static final int CACHE_SIZE = 100_000;

    @Override
    protected int maxSize() {
        return CACHE_SIZE;
    }

    @Override
    protected double compute(final int n) {
        return Gamma.digamma(n);
    }
}
