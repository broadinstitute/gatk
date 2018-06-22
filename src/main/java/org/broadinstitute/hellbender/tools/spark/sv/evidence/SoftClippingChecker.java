package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

import java.util.List;
import java.util.ListIterator;

public class SoftClippingChecker {
    public static final int MIN_SOFT_CLIP_LEN = 30; // minimum length of an interesting soft clip
    public static final byte MIN_QUALITY = 15; // minimum acceptable quality in a soft-clip window
    public static final int MAX_LOW_QUALITY_SCORES = 3; // maximum # of low quality base calls in soft-clip window

    public static boolean isClipped( final List<CigarElement> cigarElements,
                                     final byte[] quals ) {
        return hasInitialSoftClip(cigarElements, quals) || hasFinalSoftClip(cigarElements, quals);
    }

    public static boolean hasInitialSoftClip( final List<CigarElement> cigarElements, final byte[] quals ) {
        final ListIterator<CigarElement> itr = cigarElements.listIterator();
        if ( !itr.hasNext() ) return false;

        CigarElement firstEle = itr.next();
        if ( firstEle.getOperator() == CigarOperator.HARD_CLIP && itr.hasNext() ) {
            firstEle = itr.next();
        }
        final int clipStart = firstEle.getLength() - MIN_SOFT_CLIP_LEN;
        return firstEle.getOperator() == CigarOperator.SOFT_CLIP &&
                clipStart >= 0 &&
                isHighQualityRegion(quals, clipStart);
    }

    public static boolean hasFinalSoftClip( final List<CigarElement> cigarElements, final byte[] quals ) {
        final ListIterator<CigarElement> itr = cigarElements.listIterator(cigarElements.size());
        if ( !itr.hasPrevious() ) return false;

        CigarElement lastEle = itr.previous();
        if ( lastEle.getOperator() == CigarOperator.HARD_CLIP && itr.hasPrevious() ) {
            lastEle = itr.previous();
        }
        return lastEle.getOperator() == CigarOperator.SOFT_CLIP &&
                lastEle.getLength() >= MIN_SOFT_CLIP_LEN &&
                isHighQualityRegion(quals, quals.length - lastEle.getLength());
    }

    public static boolean isHighQualityRegion( final byte[] quals, int idx ) {
        if ( quals == null || quals.length == 0 ) {
            return true;
        }

        int lowQuals = 0;
        for ( final int end = idx+MIN_SOFT_CLIP_LEN; idx != end; ++idx ) {
            if ( quals[idx] < MIN_QUALITY ) {
                lowQuals += 1;
                if ( lowQuals > MAX_LOW_QUALITY_SCORES ) return false;
            }

        }
        return true;
    }
}
