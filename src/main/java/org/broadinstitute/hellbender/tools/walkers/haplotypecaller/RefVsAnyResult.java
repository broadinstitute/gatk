package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

/**
 * Holds information about a genotype call of a single sample reference vs. any non-ref event
 */
final class RefVsAnyResult {
    /**
     * The genotype likelihoods for ref/ref ref/non-ref non-ref/non-ref
     */
    private final double[] genotypeLikelihoods;

    /**
     * AD field value for ref / non-ref
     */
    private final int[] AD_Ref_Any = new int[2];


    /**
     * Creates a new ref-vs-alt result indicating the genotype likelihood vector capacity.
     * @param likelihoodCapacity the required capacity of the likelihood array, should match the possible number of
     *                           genotypes given the number of alleles (always 2), ploidy (arbitrary) less the genotyping
     *                           model non-sense genotype count if applies.
     * @throws IllegalArgumentException if {@code likelihoodCapacity} is negative.
     */
    public RefVsAnyResult(final int likelihoodCapacity) {
        if (likelihoodCapacity < 0) {
            throw new IllegalArgumentException("likelihood capacity is negative");
        }
        genotypeLikelihoods = new double[likelihoodCapacity];
    }

    /**
     * @return Get the DP (sum of AD values)
     */
    int getDP() { return AD_Ref_Any[0] + AD_Ref_Any[1]; }

    /**
     * Return the AD fields. Returns the live array.
     */
    int[] get_AD_Ref_Any(){
        return AD_Ref_Any;
    }

    double[] getGenotypeLikelihoods(){
        return genotypeLikelihoods;
    }

    void incrementRefAD(final int by){
        AD_Ref_Any[0] += by;
    }

    void incrementNonRefAD(final int by){
        AD_Ref_Any[1] += by;
    }

    /**
     * Cap the het and hom var likelihood values by the hom ref likelihood.
     */
    void capByHomRefLikelihood() {
        final int likelihoodCount = genotypeLikelihoods.length;
        for (int i = 1; i < likelihoodCount; i++) {
            genotypeLikelihoods[i] = Math.min(genotypeLikelihoods[0], genotypeLikelihoods[i]);
        }
    }

    public void increaseGenotypeLikelihood(final int idx, final double by) {
        genotypeLikelihoods[idx] += by;
    }
}
