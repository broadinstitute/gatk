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
     * mutations (SNPs); it will not produce a comprehensive set of indels. */
    EMIT_ALL_SITES
}
