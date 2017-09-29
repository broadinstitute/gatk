package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;

/**
 * Various utility functions helping calling structural variants.
 */
public final class SVVariantDiscoveryUtils {

    /**
     * @return the number of bases of two alignment regions overlap on the locally-assembled contig they originate from.
     *          Mostly useful for computing micro-homologyForwardStrandRep.
     */
    @VisibleForTesting
    public static int overlapOnContig(final AlignedAssembly.AlignmentInterval one, final AlignedAssembly.AlignmentInterval two) {
        return Math.max(0, Math.min(one.endInAssembledContig + 1, two.endInAssembledContig + 1) - Math.max(one.startInAssembledContig, two.startInAssembledContig));
    }

    /**
     * @return the total number of hard clipped bases represented in the CIGAR.
     */
    @VisibleForTesting
    public static int getTotalHardClipping(final Cigar cigar) {
        final List<CigarElement> cigarElements = cigar.getCigarElements();
        final int sz = cigarElements.size();
        if (sz <2) { // no cigar elements or only 1 element means there cannot be any hard clipping
            return 0;
        }
        return (cigarElements.get(0).getOperator() == CigarOperator.HARD_CLIP ? cigarElements.get(0).getLength() : 0) +
                (cigarElements.get(sz - 1).getOperator() == CigarOperator.HARD_CLIP ? cigarElements.get(sz - 1).getLength() : 0);
    }

    /**
     * Returns the number of clipped bases, including both soft and hard, represented in {@code cigarAlong5to3DirectionOfContig}
     * from the start or from the end
     * @param fromStart from the start of the template or not
     * @param cigar     the {@link Cigar} to be inspected
     */
    @VisibleForTesting
    public static int getNumClippedBases(final boolean fromStart, final Cigar cigar) {
        return getNumClippedBases(fromStart, cigar.getCigarElements());
    }

    /**
     * Returns the number of clipped bases, including both soft and hard, represented in {@code cigarElements}
     * from the start or from the end
     * @param fromStart     from the start of the template or not
     * @param cigarElements the ordered {@link CigarElement}'s of a cigar
     */
    @VisibleForTesting
    public static int getNumClippedBases(final boolean fromStart, final List<CigarElement> cigarElements) {

        final int sz = cigarElements.size();
        if(sz==1) return 0; // cannot be a giant clip

        final int step = fromStart ? 1 : -1;
        int result = 0;
        int j = fromStart ? 0 : sz - 1;
        CigarElement ce = cigarElements.get(j);
        while (ce.getOperator().isClipping()) {
            result += ce.getLength();
            j += step;
            if ( j < 0 || j >= sz ) break;
            ce = cigarElements.get(j);
        }
        return result;
    }

    /**
     * @return the number of hard clipped bases as indicated in the input {@code cigarElements}, either from the beginning
     * of the list ({@code fromStart==true}) or from the end of the list ({@code fromStart==false}).
     *
     * @throws IllegalArgumentException if fails check by {@link #validateCigar(List)}
     */
    @VisibleForTesting
    public static int getNumHardClippingBases(final boolean fromStart, final List<CigarElement> cigarElements) {

        validateCigar(cigarElements);

        // "H can only be present as the first and/or last operation" according to VCF spec 4.2
        final int index = fromStart ? 0 : cigarElements.size()-1;
        final CigarElement firstElement = cigarElements.get(index);
        return firstElement.getOperator()== CigarOperator.H ? firstElement.getLength() : 0;
    }

    /**
     * @return the number of soft clipped bases as indicated in the input {@code cigarElements}, either from the beginning
     * of the list ({@code fromStart==true}) or from the end of the list ({@code fromStart==false}).
     *
     * @throws IllegalArgumentException if fails check by {@link #validateCigar(List)}
     */
    @VisibleForTesting
    public static int getNumSoftClippingBases(final boolean fromStart, final List<CigarElement> cigarElements) {

        validateCigar(cigarElements);

        // because 'H' can only be the 1st/last operation according to the spec, and also
        // "S may only have H operations between them and the ends of the CIGAR string",
        // 'S' could only be the 1st operation or the 2nd operation next to the 'H' sitting at the end
        // no two 'S' operations should sit next to each other

        final int endIndex = fromStart ? 0 : cigarElements.size()-1;
        CigarElement element = cigarElements.get(endIndex);
        if (element.getOperator().isClipping()) {
            if (element.getOperator()==CigarOperator.S) {
                return element.getLength();
            } else {
                final CigarElement mayBeSoftClipping = cigarElements.get(fromStart ? 1 : cigarElements.size()-2);
                return mayBeSoftClipping.getOperator()==CigarOperator.S ? mayBeSoftClipping.getLength() : 0;
            }
        } else {
            return 0;
        }
    }

    /**
     * Checks input list of cigar operations for:
     * <ul>
     *     <li>empty input list;</li>
     *     <li>there must be at least one alignment operation in the list;</li>
     *     <li>deletion operation cannot neighbor clipping operations;</li>
     * </ul>
     * @param cigarElements
     */
    @VisibleForTesting
    public static void validateCigar(final List<CigarElement> cigarElements) {
        Utils.validateArg(!cigarElements.isEmpty(), "Cannot parse empty list cigarElements");
        Utils.validateArg(cigarElements.stream().anyMatch(ele -> ele.getOperator().isAlignment()),
                "No alignment found in the input list of cigar operations: " + cigarElements.toString());

        int idx = findIndexOfFirstNonClippingOperation(cigarElements, true);
        Utils.validateArg(idx==0 || cigarElements.get(idx).getOperator()!=CigarOperator.D,
                "Unexpected CIGAR format with deletion neighboring clipping; cigar elements are: " + cigarElements.toString());
        idx = findIndexOfFirstNonClippingOperation(cigarElements, false);
        Utils.validateArg(idx==cigarElements.size()-1 || cigarElements.get(idx).getOperator()!=CigarOperator.D,
                "Unexpected CIGAR format with deletion neighboring clipping; cigar elements are: " + cigarElements.toString());
    }

    /**
     * Returns the index of the first non-clipping operation into the input {@code cigarElements}.
     * @param cigarElements          input list of operations to be scanned through
     * @param fromStartInsteadOfEnd  either from the start of the list or from the end of the list
     */
    @VisibleForTesting
    public static int findIndexOfFirstNonClippingOperation(final List<CigarElement> cigarElements, final boolean fromStartInsteadOfEnd) {
        int idx = 0;
        final int step;
        if (fromStartInsteadOfEnd) {
            step = 1;
        } else {
            idx = cigarElements.size()-1;
            step = -1;
        }
        while(cigarElements.get(idx).getOperator().isClipping()){
            idx += step;
        }
        return idx;
    }
}
