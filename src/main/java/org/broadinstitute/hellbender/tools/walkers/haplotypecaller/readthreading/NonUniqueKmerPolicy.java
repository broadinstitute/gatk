package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

public enum NonUniqueKmerPolicy {
    NOT_ALLOWED,
    ALLOWED,
    FAIL_OVER_ONLY,
    IN_ADDITION
}
