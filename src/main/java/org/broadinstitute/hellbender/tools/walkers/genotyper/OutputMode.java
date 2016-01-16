package org.broadinstitute.hellbender.tools.walkers.genotyper;

/**
* Describes the mode of output for the caller.
*/
public enum OutputMode {
    /** produces calls only at variant sites */
    EMIT_VARIANTS_ONLY,

    /** produces calls at variant sites and confident reference sites */
    EMIT_ALL_CONFIDENT_SITES,

    /** produces calls at any callable site regardless of confidence; this argument is intended only for point
     * mutations (SNPs) in DISCOVERY mode or generally when running in GENOTYPE_GIVEN_ALLELES mode; it will by
     * no means produce a comprehensive set of indels in DISCOVERY mode */
    EMIT_ALL_SITES
}
