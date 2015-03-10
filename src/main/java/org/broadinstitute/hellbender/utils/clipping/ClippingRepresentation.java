package org.broadinstitute.hellbender.utils.clipping;

/**
 * How should we represent a clipped bases in a read?
 */
public enum ClippingRepresentation {
    /** Clipped bases are changed to Ns */
    WRITE_NS,

    /** Clipped bases are changed to have Q0 quality score */
    WRITE_Q0S,

    /** Clipped bases are change to have both an N base and a Q0 quality score */
    WRITE_NS_Q0S,

    /**
     * Change the read's cigar string to soft clip (S, see sam-spec) away the bases.
     * Note that this can only be applied to cases where the clipped bases occur
     * at the start or end of a read.
     */
    SOFTCLIP_BASES,

    /**
     * WARNING: THIS OPTION IS STILL UNDER DEVELOPMENT AND IS NOT SUPPORTED.
     *
     * Change the read's cigar string to hard clip (H, see sam-spec) away the bases.
     * Hard clipping, unlike soft clipping, actually removes bases from the read,
     * reducing the resulting file's size but introducing an irrevesible (i.e.,
     * lossy) operation.  Note that this can only be applied to cases where the clipped
     * bases occur at the start or end of a read.
     */
    HARDCLIP_BASES,

    /**
     * Turn all soft-clipped bases into matches
     */
    REVERT_SOFTCLIPPED_BASES,
}
