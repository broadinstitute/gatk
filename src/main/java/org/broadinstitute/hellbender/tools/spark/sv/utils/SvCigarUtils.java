package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.ClippingTail;

import java.util.ArrayList;
import java.util.List;

/**
 * Various utility functions helping calling structural variants.
 */
public final class SvCigarUtils {

    /**
     * Checks input list of cigar operations for:
     * <ul>
     *     <li>empty input list;</li>
     *     <li>there must be at least one alignment operation in the list;</li>
     *     <li>deletion operation cannot neighbor clipping operations;</li>
     * </ul>
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

    /**
     * Checks the input CIGAR for assumption that operator 'D' is not immediately adjacent to clipping operators.
     * Then convert the 'I' CigarElement, if it is at either end (terminal) of the input cigar, to a corresponding 'S' operator.
     * Note that we allow CIGAR of the format '10H10S10I10M', but disallows the format if after the conversion the cigar turns into a giant clip,
     * e.g. '10H10S10I10S10H' is not allowed (if allowed, it becomes a giant clip of '10H30S10H' which is non-sense).
     *
     * @return a pair of number of clipped (hard and soft, including the ones from the converted terminal 'I') bases at the front and back of the
     *         input {@code cigarAlongInput5to3Direction}.
     *
     * @throws IllegalArgumentException when the checks as described above fail.
     */
    @VisibleForTesting
    public static List<CigarElement> checkCigarAndConvertTerminalInsertionToSoftClip(final Cigar cigar) {

        if (cigar.numCigarElements()<2 ) return cigar.getCigarElements();

        final List<CigarElement> cigarElements = new ArrayList<>(cigar.getCigarElements());
        validateCigar(cigarElements);

        final List<CigarElement> convertedList = convertInsToSoftClipFromOneEnd(cigarElements, true);
        return convertInsToSoftClipFromOneEnd(convertedList, false);
    }

    /**
     * Actually convert terminal 'I' to 'S' and in case there's an 'S' comes before 'I', compactify the two neighboring 'S' operations into one.
     *
     * @return the converted and compactified list of cigar elements
     */
    @VisibleForTesting
    public static List<CigarElement> convertInsToSoftClipFromOneEnd(final List<CigarElement> cigarElements,
                                                                    final boolean fromStart) {
        final int numHardClippingBasesFromOneEnd = CigarUtils.countClippedBases(new Cigar(cigarElements), fromStart ? ClippingTail.LEFT_TAIL : ClippingTail.RIGHT_TAIL, CigarOperator.HARD_CLIP);
        final int numSoftClippingBasesFromOneEnd = CigarUtils.countClippedBases(new Cigar(cigarElements), fromStart ? ClippingTail.LEFT_TAIL : ClippingTail.RIGHT_TAIL, CigarOperator.SOFT_CLIP);

        final int indexOfFirstNonClippingOperation;
        if (numHardClippingBasesFromOneEnd==0 && numSoftClippingBasesFromOneEnd==0) { // no clipping
            indexOfFirstNonClippingOperation = fromStart ? 0 : cigarElements.size()-1;
        } else if (numHardClippingBasesFromOneEnd==0 || numSoftClippingBasesFromOneEnd==0) { // one clipping
            indexOfFirstNonClippingOperation = fromStart ? 1 : cigarElements.size()-2;
        } else {
            indexOfFirstNonClippingOperation = fromStart ? 2 : cigarElements.size()-3;
        }

        final CigarElement element = cigarElements.get(indexOfFirstNonClippingOperation);
        if (element.getOperator() == CigarOperator.I) {

            cigarElements.set(indexOfFirstNonClippingOperation, new CigarElement(element.getLength(), CigarOperator.S));

            return compactifyNeighboringSoftClippings(cigarElements);
        } else {
            return cigarElements;
        }
    }

    /**
     * Compactify two neighboring soft clippings, one of which was converted from an insertion operation.
     * @return the compactified list of operations
     * @throws IllegalArgumentException if there's un-handled edge case where two operations neighboring each other have
     *                                  the same operator (other than 'S') but for some reason was not compactified into one
     */
    @VisibleForTesting
    public static List<CigarElement> compactifyNeighboringSoftClippings(final List<CigarElement> cigarElements) {
        final List<CigarElement> result = new ArrayList<>(cigarElements.size());
        for (final CigarElement element : cigarElements) {
            final int idx = result.size()-1;
            if (result.isEmpty() || result.get(idx).getOperator()!=element.getOperator()) {
                result.add(element);
            } else {
                Utils.validateArg(result.get(idx).getOperator()==CigarOperator.S && element.getOperator()==CigarOperator.S,
                        "Seeing new edge case where two neighboring operations are having the same operator: " + cigarElements.toString());
                result.set(idx, new CigarElement(result.get(idx).getLength()+element.getLength(), CigarOperator.S));
            }
        }
        return result;
    }

}
