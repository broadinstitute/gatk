package org.broadinstitute.hellbender.tools.walkers.genotyper;

public class FRDBQDUtils {

    /**
     * These two methods are used to calculate the homopolymer base phred scaled adjustment in BQD. This code is taken from DRAGEN.
     *
     * @param paddedReference       reference to check for homopolymer span
     * @param offsetForRefIntoEvent offset of the base upon which to make a call
     */
    public static double computeForwardHomopolymerAdjustment(final byte[] paddedReference, final int offsetForRefIntoEvent, final byte errorBase) {
        int length = 1;
        while(length < 4) {
            if (errorBase != paddedReference[offsetForRefIntoEvent - length]) {
                length--;
                break;
            }
            length++;
        }
        return DRAGENGenotypesModel.BQD_HOMOPOLYMER_PHRED_ADJUSTMENT_FACTOR * length;
    }
    public static double computeReverseHomopolymerAdjustment(final byte[] paddedReference, final int offsetForRefIntoEvent, final byte errorBase) {
        int length = 1;
        while (length < 4) {
            if (errorBase != paddedReference[offsetForRefIntoEvent + length]) {
                length--;
                break;
            }
            length++;
        }
        return DRAGENGenotypesModel.BQD_HOMOPOLYMER_PHRED_ADJUSTMENT_FACTOR * length;
    }

}
