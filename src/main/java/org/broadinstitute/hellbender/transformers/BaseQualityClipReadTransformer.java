package org.broadinstitute.hellbender.transformers;

import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;

/**
 * Clips reads on both ends using base quality scores
 */
public class BaseQualityClipReadTransformer implements ReadTransformer {
    public static final long serialVersionUID = 1L;
    private int qTrimmingThreshold = 15;

    public BaseQualityClipReadTransformer(int trim_thresh) {
        qTrimmingThreshold = trim_thresh;
    }

    /**
     * Clip bases on the right end of the read from
     * <p/>
     * argmax_x{ \sum{i = x + 1}^l (qTrimmingThreshold - qual).
     * <p/>
     * Walk through the read from the right end (in machine cycle order) to the left end, calculating the
     * running sum of qTrimmingThreshold - qual.  While we do this, we track the maximum value of this
     * sum where the delta > 0.  After the loop, clipPoint is either -1 (don't do anything) or the
     * clipping index in the read (from the end).
     *
     * Repeat in reverse to clip the left end of the read.
     *
     */
    @Override
    public GATKRead apply(GATKRead read) {
        GATKRead readClippedRightEnd = clipReadRightEnd(read);
        return clipReadLeftEnd(readClippedRightEnd);
    }

    public int getRightClipPoint( final byte[] quals ) {
        int clipSum = 0, lastMax = -1, clipPoint = -1; // -1 means no clip
        final int readLength = quals.length;
        for (int i = readLength - 1; i >= 0; i--) {
            clipSum += (qTrimmingThreshold - quals[i]);
            if (clipSum >= 0 && (clipSum >= lastMax)) {
                lastMax = clipSum;
                clipPoint = i;
            }
        }
        return clipPoint;
    }

    public int getLeftClipPoint( final byte[] quals ) {
        int clipSum = 0, lastMax = -1, clipPoint = -1; // -1 means no clip
        final int readLength = quals.length;
        for (int i = 0; i < readLength; i++) {
            clipSum += (qTrimmingThreshold - quals[i]);
            if (clipSum >= 0 && (clipSum >= lastMax)) {
                lastMax = clipSum;
                clipPoint = i;
            }
        }
        return clipPoint;
    }

    private GATKRead clipReadRightEnd(GATKRead read) {
        final byte[] quals = read.getBaseQualities();
        final int clipPoint = getRightClipPoint(quals);

        if (clipPoint != -1) {
            final byte[] newBases = Arrays.copyOf(read.getBases(), clipPoint);
            final byte[] newQuals = Arrays.copyOf(quals, clipPoint);
            final GATKRead clippedRead = read.copy();
            clippedRead.setBaseQualities(newQuals);
            clippedRead.setBases(newBases);
            return clippedRead;
        } else {
            return read;
        }
    }

    private GATKRead clipReadLeftEnd(GATKRead read) {
        final byte[] quals = read.getBaseQualities();
        final int clipPoint = getLeftClipPoint(quals);

        if (clipPoint != -1) {
            final int readLength = read.getLength();
            final byte[] newBases = Arrays.copyOfRange(read.getBases(), clipPoint+1, readLength);
            final byte[] newQuals = Arrays.copyOfRange(quals, clipPoint+1, readLength);
            final GATKRead clippedRead = read.copy();
            clippedRead.setBaseQualities(newQuals);
            clippedRead.setBases(newBases);
            return clippedRead;
        } else {
            return read;
        }
    }

}


