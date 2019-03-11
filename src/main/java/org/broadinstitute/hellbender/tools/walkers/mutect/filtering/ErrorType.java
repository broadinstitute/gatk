package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

/**
 * Errors in Mutect2 fall into three major categories -- technical artifacts that depend on (usually hidden) features and do not
 * follow the independent reads assumption of the somatic likelihoods model, non-somatic variants such as germline
 * mutations and contamination, and sequencing errors that are captured by the base qualities and the somatic likelihoods model.
 *
 * It is helpful to divide errors into categories for several reasons.  First, different filters within a category are not
 * independent, and so if the probabilities of two error types within a category are p1 and p2 we may not say that the overall
 * probability is 1 - (1-p1)*(1-p2).  Instead we more conservatively take the maximum of p1 and p2.  Between categories,
 * however, we may treat error probabilities as independent.
 */
public enum ErrorType {
    ARTIFACT, NON_SOMATIC, SEQUENCING
}
