package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

public class RefVsSpanningVsAnyResult extends RefVsAnyResult {
    public int spanningDepth = 0;
    public int ployidy = 0;

    /**
     * Creates a new ref-vs-alt result indicating the genotype likelihood vector capacity.
     *
     * @param likelihoodCapacity the required capacity of the likelihood array, should match the possible number of
     *                           genotypes given the number of alleles (always 2), ploidy (arbitrary) less the genotyping
     *                           model non-sense genotype count if applies.
     * @throws IllegalArgumentException if {@code likelihoodCapacity} is negative.
     */
    public RefVsSpanningVsAnyResult(int likelihoodCapacity, int ployidy) {
        super(likelihoodCapacity);
        this.ployidy = ployidy;
    }

    /**
     * @return Get the DP (sum of AD values)
     */
    @Override
    int getDP() {
        return refDepth + spanningDepth + nonRefDepth;
    }

    /**
     * Return the AD fields. Returns a newly allocated array every time.
     */
    @Override
    int[] getAD() {
        return new int[]{refDepth, spanningDepth, nonRefDepth};
    }

    /**
     * Returns (a copy of) the array of genotype likelihoods
     * Caps the het and hom var likelihood values by the hom ref likelihood.
     * The capping is done on the fly.
     *
     * This allows any combination of ref-star alleles to be uncapped
     */
    double[] getGenotypeLikelihoodsCappedByHomRefLikelihood() {
        final double[] output = new double[genotypeLikelihoods.length];
        double minScore = 0;
        for (int i = 0; i <= ployidy; i++) {
            output[i] = genotypeLikelihoods[i];
            minScore = Math.min(minScore,genotypeLikelihoods[i]);
        }
        for (int i = ployidy + 1; i < genotypeLikelihoods.length; i++) {
            output[i] = Math.min(genotypeLikelihoods[i], minScore);
        }
        return output;
    }
}
