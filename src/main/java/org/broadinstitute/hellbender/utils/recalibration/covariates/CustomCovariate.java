package org.broadinstitute.hellbender.utils.recalibration.covariates;

/**
 * An interface to classify Covariate classes into:
 *
 * 1. Required (ReadGroup, QualityScore)
 * 2. Standard (Cycle, Context)
 * 3. Custom (any covariates defined by the user e.g. RepeatLength)
 *
 * We call 2 and 3 together the "additional" covariates.
 */
public interface CustomCovariate extends Covariate {
}
