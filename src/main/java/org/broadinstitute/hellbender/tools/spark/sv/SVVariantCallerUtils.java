package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

import java.util.List;

/**
 * Various utility functions helping calling structural variants.
 */
final class SVVariantCallerUtils {

    /**
     * @return the number of bases of two alignment regions overlap on the locally-assembled contig they originate from.
     *          Mostly useful for computing micro-homologyForwardStrandRep.
     */
    @VisibleForTesting
    static int overlapOnContig(final AlignmentRegion one, final AlignmentRegion two) {
        return Math.max(0, Math.min(one.endInAssembledContig + 1, two.endInAssembledContig + 1) - Math.max(one.startInAssembledContig, two.startInAssembledContig));
    }

    /**
     * @return the total number of hard clipped bases represented in the CIGAR.
     */
    @VisibleForTesting
    static int getTotalHardClipping(final Cigar cigar) {
        final List<CigarElement> cigarElements = cigar.getCigarElements();
        final int sz = cigarElements.size();
        if (sz <2) { // no cigar elements or only 1 element means there cannot be any hard clipping
            return 0;
        }
        return (cigarElements.get(0).getOperator() == CigarOperator.HARD_CLIP ? cigarElements.get(0).getLength() : 0) +
                (cigarElements.get(sz - 1).getOperator() == CigarOperator.HARD_CLIP ? cigarElements.get(sz - 1).getLength() : 0);
    }

    /**
     * Returns the number of clipped bases, including both soft and hard, represented in {@code forwardStrandCigar}
     * from the start or from the end
     * @param fromStart             from the start of the template or not
     * @param forwardStrandCigar    the {@link Cigar} in its forward strand representation
     */
    @VisibleForTesting
    static int getNumClippedBases(final boolean fromStart, final Cigar forwardStrandCigar) {
        final List<CigarElement> elements = forwardStrandCigar.getCigarElements();
        final int sz = elements.size();
        if(sz==1) return 0; // cannot be a giant clip

        final int step = fromStart ? 1 : -1;
        int result = 0;
        int j = fromStart ? 0 : sz - 1;
        CigarElement ce = elements.get(j);
        while (ce.getOperator().isClipping()) {
            result += ce.getLength();
            j += step;
            if ( j < 0 || j >= sz ) break;
            ce = elements.get(j);
        }
        return result;
    }
}
