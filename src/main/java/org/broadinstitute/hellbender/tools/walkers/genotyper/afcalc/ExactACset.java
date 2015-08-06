package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.Arrays;

/**
 * This class represents a column in the Exact AC calculation matrix
 */
public final class ExactACset {
    // the counts of the various alternate alleles which this column represents
    private final ExactACcounts acCounts;

    // the column of the matrix
    private final double[] log10Likelihoods;

    int sum = -1;

    public ExactACset(final int size, final ExactACcounts ACcounts) {
        this.acCounts = ACcounts;
        log10Likelihoods = new double[size];
        Arrays.fill(log10Likelihoods, Double.NEGATIVE_INFINITY);
    }

    /**
     * sum of all the non-reference alleles
     */
    public int getACsum() {
        if ( sum == -1 ) {
            sum = (int) MathUtils.sum(acCounts.getCounts());
        }
        return sum;
    }

    @Override
    public boolean equals(final Object obj) {
        return (obj instanceof ExactACset) && acCounts.equals(((ExactACset) obj).acCounts);
    }

    @Override
    public int hashCode() {
        return acCounts.hashCode();
    }

    public ExactACcounts getACcounts() {
        return acCounts;
    }

    public double[] getLog10Likelihoods() {
        return log10Likelihoods;
    }
}
