package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceResult;

/**
 * Holds information about a genotype call of a single sample reference vs. any non-ref event
 *
 * IMPORTANT PERFORMANCE NOTE!!! Allowing direct field access (within this class only) speeds up
 * the HaplotypeCaller (and thusly Mutect2) by ~10% vs. accessing the fields indirectly via setters, as seen in a profiler.
 */
final class SomaticRefVsAnyResult extends ReferenceConfidenceResult {
    /**
     * The log odds for each alt allele versus the reference
     *
     * Fields are visible because direct field access for this particular class has a major performance
     * impact on the HaplotypeCaller, as noted above, and the class itself is nested within
     * ReferenceConfidenceModel anyway.
     */
    PerAlleleCollection<Double> lods;

    /**
     * Creates a new ref-vs-alt result indicating the genotype likelihood vector capacity.
     */
    SomaticRefVsAnyResult() {
        lods = new PerAlleleCollection<>(PerAlleleCollection.Type.ALT_ONLY);
    }

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
