package org.broadinstitute.hellbender.tools.walkers.genotyper;

/**
* Describes the mode of output for the caller.
*/
public enum OutputMode {
    /** produces calls only at variant sites */
    EMIT_VARIANTS_ONLY,

    /** produces calls at variant sites and confident reference sites */
    EMIT_ALL_CONFIDENT_SITES,

    /** Produces calls at any region over the activity threshold regardless of confidence. On occasion, this will output
     * HOM_REF records where no call could be confidently made. This does not necessarily output calls for all sites in
     * a region. This argument is intended only for point mutations (SNPs); it will not produce a comprehensive set of
     * indels. */
    EMIT_ALL_ACTIVE_SITES
}
