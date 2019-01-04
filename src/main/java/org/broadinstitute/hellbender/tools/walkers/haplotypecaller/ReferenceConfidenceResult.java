package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

/**
 * Holds information about a genotype call of a single sample reference vs. any non-ref event
 */
public abstract class ReferenceConfidenceResult {
    public int refDepth = 0;
    public int nonRefDepth = 0;
    //this is abstract because it needs some kind of structure to hold likelihoods, depending on the application

    /**
     * @return Get the DP (sum of AD values)
     */
    int getDP() {
        return refDepth + nonRefDepth;
    }

    /**
     * Return the AD fields. Returns a newly allocated array every time.
     */
    int[] getAD() {
        return new int[]{refDepth, nonRefDepth};
    }
}
