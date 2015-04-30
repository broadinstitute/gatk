package org.broadinstitute.hellbender.utils.fragments;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.QualityUtil;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * Set of utilities that deal with reads that came (or are assumed to come) from the same fragment.
 */
public final class FragmentUtils {

    private FragmentUtils() {} // no instances

    public final static double DEFAULT_PCR_ERROR_RATE = 1e-4;
    public final static int DEFAULT_PCR_ERROR_QUAL = QualityUtil.getPhredScoreFromErrorProbability(DEFAULT_PCR_ERROR_RATE);
    public final static int HALF_OF_DEFAULT_PCR_ERROR_QUAL = DEFAULT_PCR_ERROR_QUAL / 2;

    /**
     * Fix two overlapping reads from the same fragment by adjusting base qualities, if possible.
     *
     * firstRead and secondRead must be part of the same fragment (though this isn't checked).  Looks
     * at the bases and alignment, and tries its best to create adjusted base qualities so that the observations
     * are not treated independently.
     *
     * This method does NOT assumes that firstRead starts before secondRead (according to their soft clipped starts).
     * Rather it checks which read starts first and calls {@link adjustQualsOfOverlappingPairedFragments} with the appropriate order.
     */
    public static void adjustQualsOfOverlappingPairedFragments(final Pair<SAMRecord, SAMRecord> overlappingPair) {
        final SAMRecord firstRead = overlappingPair.getLeft();
        final SAMRecord secondRead = overlappingPair.getRight();

        if (ReadUtils.getSoftStart(secondRead) < ReadUtils.getSoftStart(firstRead) ) {
            adjustQualsOfOverlappingPairedFragments(secondRead, firstRead);
        } else {
            adjustQualsOfOverlappingPairedFragments(firstRead, secondRead);
        }
    }

    /**
     * Fix two overlapping reads from the same fragment by adjusting base qualities, if possible
     *
     * firstRead and secondRead must be part of the same fragment (though this isn't checked).  Looks
     * at the bases and alignment, and tries its best to create adjusted base qualities so that the observations
     * are not treated independently.
     *
     * Assumes that firstRead starts before secondRead (according to their soft clipped starts)
     *
     * @param clippedFirstRead the left most read
     * @param clippedSecondRead the right most read
     */
    public static void adjustQualsOfOverlappingPairedFragments(final SAMRecord clippedFirstRead, final SAMRecord clippedSecondRead) {
        if ( clippedFirstRead == null ) {
            throw new IllegalArgumentException("clippedFirstRead cannot be null");
        }
        if ( clippedSecondRead == null ) {
            throw new IllegalArgumentException("clippedSecondRead cannot be null");
        }
        if ( ! clippedFirstRead.getReadName().equals(clippedSecondRead.getReadName()) ) {
            throw new IllegalArgumentException("attempting to merge two reads with different names " + clippedFirstRead + " and " + clippedSecondRead);
        }

        // don't adjust fragments that do not overlap
        final boolean overlap = clippedFirstRead.getReferenceIndex().equals(clippedSecondRead.getReferenceIndex()) && clippedFirstRead.getAlignmentEnd() >= clippedSecondRead.getAlignmentStart();
        if ( !overlap ) {
            return;
        }

        final Pair<Integer, Boolean> pair = ReadUtils.getReadCoordinateForReferenceCoordinate(clippedFirstRead, clippedSecondRead.getAlignmentStart());
        final int firstReadStop = ( pair.getRight() ? pair.getLeft() + 1 : pair.getLeft() );
        final int numOverlappingBases = Math.min(clippedFirstRead.getReadLength() - firstReadStop, clippedSecondRead.getReadLength());

        final byte[] firstReadBases = clippedFirstRead.getReadBases();
        final byte[] firstReadQuals = clippedFirstRead.getBaseQualities();
        final byte[] secondReadBases = clippedSecondRead.getReadBases();
        final byte[] secondReadQuals = clippedSecondRead.getBaseQualities();

        for ( int i = 0; i < numOverlappingBases; i++ ) {
            final int firstReadIndex = firstReadStop + i;
            final byte firstReadBase = firstReadBases[firstReadIndex];
            final byte secondReadBase = secondReadBases[i];

            if ( firstReadBase == secondReadBase ) {
                firstReadQuals[firstReadIndex] = (byte) Math.min(firstReadQuals[firstReadIndex], HALF_OF_DEFAULT_PCR_ERROR_QUAL);
                secondReadQuals[i] = (byte) Math.min(secondReadQuals[i], HALF_OF_DEFAULT_PCR_ERROR_QUAL);
            } else {
                // TODO -- use the proper statistical treatment of the quals from DiploidSNPGenotypeLikelihoods.java
                firstReadQuals[firstReadIndex] = 0;
                secondReadQuals[i] = 0;
            }
        }

        clippedFirstRead.setBaseQualities(firstReadQuals);
        clippedSecondRead.setBaseQualities(secondReadQuals);
    }
}
