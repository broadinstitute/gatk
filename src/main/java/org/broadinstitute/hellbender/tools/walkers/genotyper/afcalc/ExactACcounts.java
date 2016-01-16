package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import java.util.Arrays;

/**
 * Stores a vector of counts. It's a thin wrapper around int[] to give it a cached hashcode and equals.
 */
public final class ExactACcounts {
    private final int[] counts;
    private int hashcode = -1;

    /**
     * Note: this constructor does not make a copy of the argument and stores a live pointer to the given array.
     * Callers must make sure the array is not mutated to maintain semantics of hashcode.
     * The array must be not null and longer than 0 elements.
     */
    public ExactACcounts(final int[] counts) {
        if (counts == null || counts.length == 0){
            throw new IllegalArgumentException("counts should not be null or empty");
        }
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
            hashcode = Arrays.hashCode(counts);
        }
        return hashcode;
    }

    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        sb.append(counts[0]);
        for ( int i = 1; i < counts.length; i++ ) {
            sb.append("/");
            sb.append(counts[i]);
        }
        return sb.toString();
    }
}
