package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.stream.DoubleStream;

/**
 * Holds information about a genotype call of a single sample reference vs. any non-ref event
 *
 * IMPORTANT PERFORMANCE NOTE!!! Allowing direct field access (within this class only) speeds up
 * the HaplotypeCaller by ~10% vs. accessing the fields indirectly via setters, as seen in a profiler.
 */
public final class RefVsAnyResult extends ReferenceConfidenceResult {
    /**
     * The genotype likelihoods for ref/ref ref/non-ref non-ref/non-ref
     *
     * Fields are visible because direct field access for this particular class has a major performance
     * impact on the HaplotypeCaller, as noted above, and the class itself is nested within
     * ReferenceConfidenceModel anyway.
     */
    final double[] genotypeLikelihoods;

    int[] finalPhredScaledGenotypeLikelihoods;

    /**
     * Creates a new ref-vs-alt result indicating the genotype likelihood vector capacity.
     * @param likelihoodCapacity the required capacity of the likelihood array, should match the possible number of
     *                           genotypes given the number of alleles (always 2), ploidy (arbitrary) less the genotyping
     *                           model non-sense genotype count if applies.
     * @throws IllegalArgumentException if {@code likelihoodCapacity} is negative.
     */
    public RefVsAnyResult(final int likelihoodCapacity) {
        ParamUtils.isPositiveOrZero(likelihoodCapacity, "likelihood capacity is negative");
        genotypeLikelihoods = new double[likelihoodCapacity];
        finalPhredScaledGenotypeLikelihoods = new int[likelihoodCapacity];
    }

    /**
     * Returns (a copy of) the array of genotype likelihoods
     * Caps the het and hom var likelihood values by the hom ref likelihood.
     * The capping is done on the fly.
     */
    double[] getGenotypeLikelihoodsCappedByHomRefLikelihood() {
        return DoubleStream.of(genotypeLikelihoods).map(d -> Math.min(d, genotypeLikelihoods[0])).toArray();
    }
}
