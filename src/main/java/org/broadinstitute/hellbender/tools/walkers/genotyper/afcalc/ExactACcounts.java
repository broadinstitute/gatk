package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import java.util.Arrays;

/**
 * Stores a vector of counts. It's a thin wrapper around int[] to give it a hash and equals.
 */
public final class ExactACcounts {
    private final int[] counts;
    private int hashcode = -1;

    /**
     * Note: this constructor does not make a copy of the argument and stores a live pointer to the given array.
     * Callers must make sure the array is not mutated to maintain semantics of hashcode.
     */
    public ExactACcounts(final int[] counts) {
        this.counts = counts;
    }

    public int[] getCounts() {
        return counts;
    }

    @Override
    public boolean equals(final Object obj) {
        return (obj instanceof ExactACcounts) && Arrays.equals(counts, ((ExactACcounts) obj).counts);
    }

    @Override
    public int hashCode() {
        if ( hashcode == -1 ) {
            hashcode = Arrays.hashCode(getCounts());
        }
        return hashcode;
    }

    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        sb.append(getCounts()[0]);
        for ( int i = 1; i < getCounts().length; i++ ) {
            sb.append("/");
            sb.append(getCounts()[i]);
        }
        return sb.toString();
    }
}
