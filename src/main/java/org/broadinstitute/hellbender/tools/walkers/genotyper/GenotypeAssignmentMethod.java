package org.broadinstitute.hellbender.tools.walkers.genotyper;

/**
 * Created by davidben on 6/10/16.
 */
public enum GenotypeAssignmentMethod {
    /**
     * set all of the genotype GT values to NO_CALL
     */
    SET_TO_NO_CALL,

    /**
     * Use the subsetted PLs to greedily assigned genotypes
     */
    USE_PLS_TO_ASSIGN,

    /**
     * do not even bother changing the GTs
     */
    DO_NOT_ASSIGN_GENOTYPES
}
