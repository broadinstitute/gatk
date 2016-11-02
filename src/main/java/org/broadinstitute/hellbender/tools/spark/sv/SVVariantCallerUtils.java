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

    @VisibleForTesting
    static int overlapOnContig(final AlignmentRegion one, final AlignmentRegion two) {
        return Math.max(0, Math.min(one.endInAssembledContig + 1, two.endInAssembledContig + 1) - Math.max(one.startInAssembledContig, two.startInAssembledContig));
    }

    /**
     * Test if a {@link BreakpointAllele} is an {@link org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.BreakpointAlleleInversion}
     */
    @VisibleForTesting
    static boolean isInversion(final BreakpointAllele allele) {
        try {
            final BreakpointAllele.BreakpointAlleleInversion invAllele = (BreakpointAllele.BreakpointAlleleInversion) allele;
            return invAllele.leftAlignedLeftBreakpoint.getContig().equals(invAllele.leftAlignedRightBreakpoint.getContig())
                    &&
                    (invAllele.getInversionType() == BreakpointAllele.BreakpointAlleleInversion.InversionType.INV_3_TO_5 || invAllele.getInversionType() == BreakpointAllele.BreakpointAlleleInversion.InversionType.INV_5_TO_3);
        } catch (final ClassCastException ccex) {
            return false;
        }
    }

    @VisibleForTesting
    static int getTotalHardClipping(final Cigar cigar) {
        final List<CigarElement> cigarElements = cigar.getCigarElements();
        final int sz = cigarElements.size();
        if (sz == 0) {
            return 0;
        }
        return (cigarElements.get(0).getOperator() == CigarOperator.HARD_CLIP ? cigarElements.get(0).getLength() : 0) +
                (sz == 1 || cigarElements.get(sz - 1).getOperator() == CigarOperator.HARD_CLIP ? cigarElements.get(sz - 1).getLength() : 0);
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
