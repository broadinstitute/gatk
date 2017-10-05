package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * Various utility functions helping calling structural variants.
 */
public final class SvCigarUtils {

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

    @VisibleForTesting
    public static int getUnclippedReadLength(final Cigar cigar) {
        validateCigar(cigar.getCigarElements());
        return cigar.getCigarElements().stream()
                .mapToInt(element ->
                        element.getOperator().isClipping() || element.getOperator().consumesReadBases() ? element.getLength() : 0
                )
                .sum();
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
        final int numHardClippingBasesFromOneEnd = getNumHardClippingBases(fromStart, cigarElements);
        final int numSoftClippingBasesFromOneEnd = getNumSoftClippingBases(fromStart, cigarElements);

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

    /**
     * Computes the corresponding distance needs to be walked on the reference, given the Cigar and distance walked on the read.
     * @param cigarAlong5To3DirOfRead   cigar along the 5-3 direction of read (when read is mapped to reverse strand, bwa mem output cigar should be inverted)
     * @param startInclusiveOnRead      start position (1-based) on the read (note it should not count the hard clipped bases, as usual)
     * @param distanceOnRead            distance to walk on the read
     * @return                          corresponding walk distance on reference
     * @throws IllegalArgumentException if input cigar contains padding operation or 'N', or
     *                                  either of startInclusive or distance is non-positive, or
     *                                  startInclusive + distance -1 is longer than the read as suggested by the cigar
     */
    @VisibleForTesting
    public static int computeAssociatedDistOnRef(final Cigar cigarAlong5To3DirOfRead, final int startInclusiveOnRead,
                                                 final int distanceOnRead) {

        final int endInclusive = startInclusiveOnRead + distanceOnRead - 1;
        Utils.validateArg(startInclusiveOnRead>0 && distanceOnRead > 0,
                "start position (" + startInclusiveOnRead + ") or distance (" + distanceOnRead + ") is non-positive.");
        final List<CigarElement> cigarElements = cigarAlong5To3DirOfRead.getCigarElements();
        Utils.validateArg(cigarElements.stream().noneMatch(ce -> ce.getOperator().isPadding() || ce.getOperator().equals(CigarOperator.N)),
                "cigar contains padding, which is currently unsupported; cigar: " + TextCigarCodec.encode(cigarAlong5To3DirOfRead));
        Utils.validateArg(cigarElements.stream()
                        .mapToInt(ce -> ce.getOperator().consumesReadBases() ? ce.getLength() : 0).sum() >= endInclusive,
                "start location (" + startInclusiveOnRead + ") and walking distance (" + distanceOnRead +
                        ") would walk out of the read, indicated by cigar " + TextCigarCodec.encode(cigarAlong5To3DirOfRead));

        int readBasesConsumed = 0;

        // skip first several elements that give accumulated readBasesConsumed below startInclusiveOnRead
        int idx = 0;
        CigarElement currEle = cigarElements.get(idx);
        while (readBasesConsumed + (currEle.getOperator().consumesReadBases() ? currEle.getLength() : 0) < startInclusiveOnRead) {
            readBasesConsumed += currEle.getOperator().consumesReadBases() ? currEle.getLength() : 0;
            currEle = cigarElements.get(++idx);
        }
        int refWalkDist = 0;
        int readWalked = 0;
        while (idx != cigarElements.size()) {
            currEle = cigarElements.get(idx);
            final int skip = Math.max(0, startInclusiveOnRead - readBasesConsumed - 1);
            final int effectiveLen = currEle.getLength() - skip;

            if (currEle.getOperator().consumesReadBases()) {
                if (readWalked + effectiveLen < distanceOnRead) { // hasn't walked enough yet on read
                    readWalked += effectiveLen;
                    refWalkDist += currEle.getOperator().consumesReferenceBases() ? effectiveLen : 0;
                    readBasesConsumed += currEle.getOperator().consumesReferenceBases() ? currEle.getLength() : 0;
                } else { // would be walked enough on read
                    refWalkDist += currEle.getOperator().consumesReferenceBases() ? distanceOnRead - readWalked : 0;
                    readWalked = distanceOnRead;
                    break;
                }
            } else {
                refWalkDist += currEle.getOperator().consumesReferenceBases() ? effectiveLen : 0;
                readBasesConsumed += currEle.getOperator().consumesReferenceBases() ? currEle.getLength() : 0;
            }

            ++idx;
        }

        return refWalkDist;
    }

    /**
     * Computes the corresponding distance needs to be walked on the read, given the Cigar and distance walked on the reference.
     * @param cigarAlong5To3DirOfRead   cigar along the 5-3 direction of read (when read is mapped to reverse strand, bwa mem output cigar should be inverted)
     * @param startInclusiveOnRead      start position (1-based) on the read (note it should not count the hard clipped bases, as usual)
     * @param distOnRef                 distance to walk on the reference
     * @param walkBackward              whether to walk backwards along the read or not
     * @return                          corresponding walk distance on read (always positive)
     * @throws IllegalArgumentException if input cigar contains padding operation or 'N', or
     *                                  either of {@code startInclusiveOnRead} or distance is non-positive, or
     *                                  {@code startInclusiveOnRead} is larger than read length, or
     *                                  requested reference walk distance is longer than the total read bases in cigar, or
     *                                  computed read walk distance would "walk off" the read
     */
    @VisibleForTesting
    public static int computeAssociatedDistOnRead(final Cigar cigarAlong5To3DirOfRead, final int startInclusiveOnRead,
                                                  final int distOnRef, final boolean walkBackward) {

        Utils.validateArg(distOnRef > 0 && startInclusiveOnRead > 0,
                "start position (" + startInclusiveOnRead + ") or distance (" + distOnRef + ") is non-positive.");

        final List<CigarElement> cigarElementsInOrderOfWalkingDir = (walkBackward ? CigarUtils.invertCigar(cigarAlong5To3DirOfRead): cigarAlong5To3DirOfRead).getCigarElements();
        Utils.validateArg(cigarElementsInOrderOfWalkingDir.stream().noneMatch(ce -> ce.getOperator().isPadding() || ce.getOperator().equals(CigarOperator.N)),
                "cigar contains padding, which is currently unsupported; cigar: " + TextCigarCodec.encode(cigarAlong5To3DirOfRead));
        final int readUnclippedLength = getUnclippedReadLength(cigarAlong5To3DirOfRead);
        Utils.validateArg(readUnclippedLength >= startInclusiveOnRead,
                "given start location on read (" + startInclusiveOnRead + ") is higher than read unclipped length (" + readUnclippedLength+ "), cigar: " + TextCigarCodec.encode(cigarAlong5To3DirOfRead));
        final int totalRefLen = cigarElementsInOrderOfWalkingDir.stream().mapToInt(ce -> ce.getOperator().consumesReferenceBases() ? ce.getLength() : 0).sum();
        Utils.validateArg(totalRefLen >= distOnRef,
                "given walking distance on reference (" + distOnRef + ") would be longer than the total number (" +
                        + totalRefLen + ") of reference bases spanned by the cigar, indicated by cigar " + TextCigarCodec.encode(cigarAlong5To3DirOfRead));

        final int readLength = cigarElementsInOrderOfWalkingDir.stream().mapToInt(ce -> ce.getOperator().consumesReadBases() ? ce.getLength() : 0).sum();
        final int effectiveReadStartInclusive = walkBackward ? readLength - startInclusiveOnRead + 1 : startInclusiveOnRead;
        // skip first several elements that give accumulated readBasesConsumed below startInclusiveOnRead
        int idx = 0;
        int readBasesConsumed = 0;
        CigarElement currEle = cigarElementsInOrderOfWalkingDir.get(idx);
        while (readBasesConsumed + (currEle.getOperator().consumesReadBases() ? currEle.getLength() : 0) < effectiveReadStartInclusive) {
            readBasesConsumed += currEle.getOperator().consumesReadBases() ? currEle.getLength() : 0;
            currEle = cigarElementsInOrderOfWalkingDir.get(++idx);
        }

        // when we reach here, we have skipped just enough read bases to start counting ref bases or, currEle would lead us to such state
        int readWalkDist = 0;
        int refWalked = 0;
        while (idx != cigarElementsInOrderOfWalkingDir.size()) {
            currEle = cigarElementsInOrderOfWalkingDir.get(idx);
            final int skip = Math.max(0, effectiveReadStartInclusive - readBasesConsumed - 1);
            final int effectiveLen = currEle.getLength() - skip;

            if (currEle.getOperator().consumesReferenceBases()) {
                if (refWalked + effectiveLen < distOnRef) { // hasn't walked enough yet on reference
                    refWalked += effectiveLen;
                    readWalkDist += currEle.getOperator().consumesReadBases() ? effectiveLen : 0;
                    readBasesConsumed += currEle.getOperator().consumesReadBases() ? currEle.getLength() : 0;
                } else { // would be walked enough on reference
                    readWalkDist += currEle.getOperator().consumesReadBases() ? distOnRef - refWalked : 0;
                    refWalked = distOnRef;
                    break;
                }
            } else {
                readWalkDist += currEle.getOperator().consumesReadBases() ? effectiveLen : 0;
                readBasesConsumed += currEle.getOperator().consumesReadBases() ? currEle.getLength() : 0;
            }
            ++idx;
        }

        if (refWalked < distOnRef)
            throw new IllegalArgumentException("Computed walk distance (start: " + startInclusiveOnRead + ", distOnRef: " +
                    distOnRef + ") on read beyond read length (" + readLength +") with cigar " + TextCigarCodec.encode(cigarAlong5To3DirOfRead));

        return readWalkDist;
    }

}
