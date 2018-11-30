package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;

/**
 * Classes that implement this interface have the ability to evaluate how likely it is that a site is "active"
 * (contains potential real variation).
 */
@FunctionalInterface
public interface AssemblyRegionEvaluator {

    /**
     * Given a pileup over a single locus, returns an ActivityProfileState containing the probability (0.0 to 1.0) that
     * the locus is an "active" site.
     *
     * @param locusPileup reads pileup to examine
     * @param referenceContext reference base overlapping the pileup locus
     * @param featureContext features overlapping the pileup locus
     * @return ActivityProfileState containing the probability between 0.0 and 1.0 that the site is active
     */
    ActivityProfileState isActive( final AlignmentContext locusPileup, final ReferenceContext referenceContext, final FeatureContext featureContext );

    // After trimming to fit the assembly window, throw away read stubs shorter than this length
    // if we don't, the several bases left of reads that end just within the assembly window can
    // get realigned incorrectly.  See https://github.com/broadinstitute/gatk/issues/5060
    int MINIMUM_READ_LENGTH_AFTER_TRIMMING = 10;
}
