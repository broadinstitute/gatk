package org.broadinstitute.hellbender.utils;

public final class Log10Cache extends IntToDoubleFunctionCache {
    @Override
    protected int maxSize() {
        return Integer.MAX_VALUE;
    }

    @Override
    protected double compute(final int n) {
        return Math.log10(n);
    }
}
