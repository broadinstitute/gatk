package org.broadinstitute.hellbender.utils.fragments;

import htsjdk.samtools.util.QualityUtil;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.List;

public final class FragmentUtils {
    private FragmentUtils() {}

    public final static double DEFAULT_PCR_ERROR_RATE = 1e-4;
    public final static int DEFAULT_PCR_ERROR_QUAL = QualityUtil.getPhredScoreFromErrorProbability(DEFAULT_PCR_ERROR_RATE);
    public final static int HALF_OF_DEFAULT_PCR_ERROR_QUAL = DEFAULT_PCR_ERROR_QUAL / 2;

    /**
     * Fix two overlapping reads from the same fragment by adjusting base qualities, if possible
     *
     * firstRead and secondRead must be part of the same fragment (though this isn't checked).  Looks
     * at the bases and alignment, and tries its best to create adjusted base qualities so that the observations
     * are not treated independently.
     *
     * Assumes that firstRead starts before secondRead (according to their soft clipped starts)
     * Sets the qualities of clippedFirstRead and clippedSecondRead to mimic a merged read or
     * nothing if the algorithm cannot create a meaningful one
     *
     * @param clippedFirstRead the left most read
     * @param clippedSecondRead the right most read
     *
     */
    public static void adjustQualsOfOverlappingPairedFragments(final GATKRead clippedFirstRead, final GATKRead clippedSecondRead) {
        Utils.nonNull(clippedFirstRead);
        Utils.nonNull(clippedSecondRead);
        Utils.validateArg(clippedFirstRead.getName().equals(clippedSecondRead.getName()), () ->
                "attempting to merge two reads with different names " + clippedFirstRead + " and " + clippedSecondRead);

        // don't adjust fragments that do not overlap
        if ( clippedFirstRead.getEnd() < clippedSecondRead.getStart() || !clippedFirstRead.getContig().equals(clippedSecondRead.getContig()) ) {
            return;
        }

        final Pair<Integer, Boolean> pair = ReadUtils.getReadCoordinateForReferenceCoordinate(clippedFirstRead, clippedSecondRead.getStart());
        final int firstReadStop = ( pair.getRight() ? pair.getLeft() + 1 : pair.getLeft() );
        final int numOverlappingBases = Math.min(clippedFirstRead.getLength() - firstReadStop, clippedSecondRead.getLength());

        final byte[] firstReadBases = clippedFirstRead.getBases();
        final byte[] firstReadQuals = clippedFirstRead.getBaseQualities();
        final byte[] secondReadBases = clippedSecondRead.getBases();
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

    public static void adjustQualsOfOverlappingPairedFragments( final List<GATKRead> overlappingPair ) {
        Utils.validateArg( overlappingPair.size() == 2, () -> "Found overlapping pair with " + overlappingPair.size() + " reads, but expecting exactly 2.");

        final GATKRead firstRead = overlappingPair.get(0);
        final GATKRead secondRead = overlappingPair.get(1);

        if ( secondRead.getSoftStart() < firstRead.getSoftStart() ) {
            adjustQualsOfOverlappingPairedFragments(secondRead, firstRead);
        }
        else {
            adjustQualsOfOverlappingPairedFragments(firstRead, secondRead);
        }
    }
}
