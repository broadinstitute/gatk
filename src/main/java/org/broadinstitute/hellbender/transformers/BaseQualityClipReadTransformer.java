package org.broadinstitute.hellbender.transformers;

import org.broadinstitute.hellbender.utils.clipping.ClippingOp;
import org.broadinstitute.hellbender.utils.clipping.ClippingRepresentation;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Clips reads on both ends using base quality scores
 */
public final class BaseQualityClipReadTransformer implements ReadTransformer {
    private static final long serialVersionUID = 1L;
    private int qTrimmingThreshold = 15;

    public BaseQualityClipReadTransformer(final int trim_thresh) {
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
    public GATKRead apply(final GATKRead read) {
        final GATKRead readClippedRightEnd = clipReadRightEnd(read);
        return clipReadLeftEnd(readClippedRightEnd);
    }

    /**
     * Returns right clip point or -1 if no clip
     */
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

    /**
     * Returns left clip point or -1 if no clip
     */
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
            final ReadClipper readClipper = new ReadClipper(read);
            readClipper.addOp(new ClippingOp(clipPoint, read.getLength()));
            return readClipper.clipRead(ClippingRepresentation.HARDCLIP_BASES);
        } else {
            return read;
        }
    }

    private GATKRead clipReadLeftEnd(GATKRead read) {
        final byte[] quals = read.getBaseQualities();
        final int clipPoint = getLeftClipPoint(quals);
        if (clipPoint != -1) {
            final ReadClipper readClipper = new ReadClipper(read);
            readClipper.addOp(new ClippingOp(0, clipPoint));
            return readClipper.clipRead(ClippingRepresentation.HARDCLIP_BASES);
        } else {
            return read;
        }
    }

}


