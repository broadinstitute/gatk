package org.broadinstitute.hellbender.utils.recalibration.covariates;

/**
 * An interface to classify Covariate classes into:
 * 1. Required (ReadGroup, QualityScore)
 * 2. Standard (Cycle, Context)
 * 3. Custom (e.g. RepeatLength)
 *
 * 2 and 3 together are called "Additional" covariates in parts of the code.
 *
 */
public interface CustomCovariate extends Covariate {
}
