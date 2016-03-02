package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

/**
 * How to resolve the haplotype graph when haplotypes where generated from a mixture of different kmerSizes.
 */
public enum HeterogeneousKmerSizeResolution {

    /**
     * Combine haplotypes using a haplotype graph with the largest kmerSize amongst the ones that generated some haplotype.
     */
    COMBO_MAX,

    /**
     * Combine haplotypes using a haplotype graph with the largest kmerSize amongst the ones that generated some haplotype.
     */
    COMBO_MIN,

    /**
     * Take just the haplotypes from largest kmersize that generated any.
     */
    MAX_ONLY,

    /**
     * Take just the haplotypes from the smallest kmerSize that generated any.
     */
    @SuppressWarnings("unused")
    MIN_ONLY;

    /**
     * Indicates whether we should use the maximum kmerSize for the haplotypeGraph or not.
     *
     * @return true if we need to use the maximum, false otherwise.
     */
    public boolean useMaximum() {
        switch (this) {
            case COMBO_MAX: return true;
            case MAX_ONLY: return true;
            default: return false;
        }
    }

    /**
     * Indicates whether we should use the minimum kmerSize for the haplotypeGraph or not.
     *
     * @return true if we need to use the minimum, false otherwise.
     */
    @SuppressWarnings("unused")
    public boolean useMinimum() {
        return ! useMaximum();
    }

    /**
     * Tell whether this policy combines kmer-sizes or not.
     * @return true iff it does.
     */
    public boolean combinesKmerSizes() {
        switch (this) {
            case COMBO_MAX: return true;
            case COMBO_MIN: return true;
            default: return false;
        }

    }
}

