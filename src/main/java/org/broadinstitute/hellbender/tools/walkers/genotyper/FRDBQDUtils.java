package org.broadinstitute.hellbender.tools.walkers.genotyper;

/**
 * Collection of methods for calculating FRD and BQD corrected genotypes.
 */
final class FRDBQDUtils {

    /**
     * Calculates the homopolymer base phred scaled adjustment in BQD that is applicable to bases downstream from a given offset on the reference.
     * <p>This code has been taken and from DRAGEN with some modifications</p>
     *
     * @param paddedReference       reference to check for homopolymer span
     * @param offsetForRefIntoEvent offset of the base upon which to make a call
     */
    static double computeForwardHomopolymerAdjustment(final byte[] paddedReference, final int offsetForRefIntoEvent, final byte errorBase) {
        int length, offset;
        for (length = 0, offset = offsetForRefIntoEvent - 1; offset >= 0 && length < 4; length++, offset--) {
            if (errorBase != paddedReference[offset]) {
                break;
            }
        }
        return DRAGENGenotypesModel.BQD_HOMOPOLYMER_PHRED_ADJUSTMENT_FACTOR * length;
    }

    /**
     * Calculates the homopolymer base phred scaled adjustment in BQD that is applicable to bases upstream from a given offset on the reference.
     * <p>This code has been taken and from DRAGEN with some modifications</p>
     *
     * @param paddedReference       reference to check for homopolymer span
     * @param offsetForRefIntoEvent offset of the base upon which to make a call
     */
    static double computeReverseHomopolymerAdjustment(final byte[] paddedReference, final int offsetForRefIntoEvent, final byte errorBase) {
        int length, offset;
        for (length = 0, offset = offsetForRefIntoEvent + 1; offset < paddedReference.length && length < 4; length++, offset++) {
            if (errorBase != paddedReference[offset]) {
                break;
            }
        }
        return DRAGENGenotypesModel.BQD_HOMOPOLYMER_PHRED_ADJUSTMENT_FACTOR * length;
    }
}
