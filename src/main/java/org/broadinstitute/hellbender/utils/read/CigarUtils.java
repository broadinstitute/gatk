package org.broadinstitute.hellbender.utils.read;

import com.google.common.collect.Lists;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.utils.Tail;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public final class CigarUtils {

    private static final String SW_PAD = "NNNNNNNNNN";

    private CigarUtils(){}

    /**
     * Inverts the order of the operators in the cigar.
     * Eg 10M1D20M -> 20M1D10M
     */
    public static Cigar invertCigar (final Cigar cigar) {
        Utils.nonNull(cigar);
        final List<CigarElement>  els = new ArrayList<>(cigar.getCigarElements());
        Collections.reverse(els);
        return new Cigar(els);
    }

    /**
     * Compute the number of reference bases between the start (inclusive) and end (exclusive) cigar elements.
     * Reference bases are counted as the number of positions w.r.t. the reference spanned by an alignment to that reference.
     * For example original position = 10. cigar: 2M3I2D1M. If you remove the 2M the new starting position is 12.
     * If you remove the 2M3I it is still 12. If you remove 2M3I2D (not reasonable cigar), you will get position 14.
     */
    public static int countRefBasesAndClips(final List<CigarElement> elems, final int cigarStartIndex, final int cigarEndIndex){
        return countRefBasesAndMaybeAlsoClips(elems, cigarStartIndex, cigarEndIndex, true, true);
    }

    public static int countRefBasesAndSoftClips(final List<CigarElement> elems, final int cigarStartIndex, final int cigarEndIndex){
        return countRefBasesAndMaybeAlsoClips(elems, cigarStartIndex, cigarEndIndex, true, false);
    }

    private static int countRefBasesAndMaybeAlsoClips(final List<CigarElement> elems, final int cigarStartIndex, final int cigarEndIndex, final boolean includeSoftClips, final boolean includeHardClips) {
        Utils.nonNull(elems);
        Utils.validateArg(cigarStartIndex >= 0 && cigarEndIndex <= elems.size() && cigarStartIndex <= cigarEndIndex, () -> "invalid index:" + 0 + " -" + elems.size());

        int result = 0;
        for (final CigarElement elem : elems.subList(cigarStartIndex, cigarEndIndex)) {
            final CigarOperator op = elem.getOperator();
            if (op.consumesReferenceBases() || (includeSoftClips && op == CigarOperator.SOFT_CLIP) || (includeHardClips && op == CigarOperator.HARD_CLIP)) {
                result += elem.getLength();
            }
        }

        return result;
    }

    /**
     * Removes all clipping and padding operators from the cigar.
     */
    public static Cigar removeClipsAndPadding(final Cigar cigar) {
        Utils.nonNull(cigar, "cigar is null");
        final List<CigarElement> elements = new ArrayList<>(cigar.numCigarElements());
        for ( final CigarElement ce : cigar.getCigarElements() ) {
            if ( !isClipOrPaddingOperator(ce.getOperator()) ) {
                elements.add(ce);
            }
        }
        return new Cigar(elements);
    }

    private static boolean isClipOrPaddingOperator(final CigarOperator op) {
        return op == CigarOperator.S || op == CigarOperator.H || op == CigarOperator.P;
    }

    /**
     * Returns whether the list has any N operators.
     */
    public static boolean containsNOperator(final List<CigarElement> cigarElements) {
        return Utils.nonNull(cigarElements).stream().anyMatch(el -> el.getOperator() == CigarOperator.N);
    }

    /**
     * A good Cigar object obeys the following rules:
     *  - is valid as per SAM spec {@link Cigar#isValid(String, long)}.
     *  - has no consecutive I/D elements
     *  - does not start or end with deletions (with or without preceding clips).
     */
    public static boolean isGood(final Cigar c) {
        Utils.nonNull(c, "cigar is null");

        //Note: the string comes from the SAMRecord so it must be a wellformed CIGAR (that is, in "\*|([0-9]+[MIDNSHPX=])+" as per SAM spec).
        //We don't have to check that
        if (c.isValid(null, -1) != null){  //if it's invalid, then it's not good
            return false;
        }
        final List<CigarElement> elems = c.getCigarElements();
        return !(hasConsecutiveIndels(elems) || startsOrEndsWithDeletionIgnoringClips(elems));
    }

    /**
     * Checks if cigar has consecutive I/D elements.
     */
    private static boolean hasConsecutiveIndels(final List<CigarElement> elems) {
        boolean prevIndel = false;
        for (final CigarElement elem : elems) {
            final CigarOperator op = elem.getOperator();
            final boolean isIndel = (op == CigarOperator.INSERTION || op == CigarOperator.DELETION);
            if (prevIndel && isIndel) {
                return true;
            }
            prevIndel = isIndel;
        }
        return false;
    }

    /**
     * Checks if cigar starts with a deletion (ignoring any clips at the beginning).
     */
    private static boolean startsOrEndsWithDeletionIgnoringClips(final List<CigarElement> elems) {

        for (final boolean leftSide : new boolean[] {true, false}) {
            for (final CigarElement elem : leftSide ? elems : Lists.reverse(elems)) {
                final CigarOperator op = elem.getOperator();
                if (op == CigarOperator.DELETION) { //first non-clipping is deletion
                    return true;
                } else if (!op.isClipping()) {  // first non-clipping is non deletion
                    break;
                }
            }
        }

        return false;
    }

    /**
     * Calculate the cigar elements for this path against the reference sequence.
     *
     * This assumes that the reference and alt sequence are haplotypes derived from a de Bruijn graph or SeqGraph and have the same
     * ref source and ref sink vertices.  That is, the alt sequence start and end are assumed anchored to the reference start and end, which
     * occur at the ends of the padded assembly region.  Hence, unlike read alignment, there is no concept of a start or end coordinate here.
     * Furthermore, it is important to note that in the rare case that the alt cigar begins or ends with a deletion, we must keep the leading
     * or trailing deletion in order to maintain the original reference span of the alt haplotype.  This can occur, for example, when the ref
     * haplotype starts with N repeats of a long sequence and the alt haplotype starts with N-1 repeats.
     *
     * @param refSeq the reference sequence that all of the bases in this path should align to
     * @return a Cigar mapping this path to refSeq, or null if no reasonable alignment could be found
     */
    public static Cigar calculateCigar(final byte[] refSeq, final byte[] altSeq, final SmithWatermanAligner aligner, final SWParameters haplotypeToReferenceSWParameters, final SWOverhangStrategy strategy) {
        Utils.nonNull(refSeq, "refSeq");
        Utils.nonNull(altSeq, "altSeq");
        if ( altSeq.length == 0 ) {
            // horrible edge case from the unit tests, where this path has no bases
            return new Cigar(Collections.singletonList(new CigarElement(refSeq.length, CigarOperator.D)));
        }

        //Note: this is a performance optimization.
        // If two strings are equal (a O(n) check) then it's trivial to get CIGAR for them.
        //TODO: in addressing http://github.com/broadinstitute/gatk/issues/6863,
        // an optimization for strings with <= 2 mismatches (introduced in https://github.com/broadinstitute/gatk/pull/5466)
        // was reverted, as it assumes SW parameters that prefer 2 substitutions over 1 insertion + 1 deletion;
        // we could restore the optimization if appropriate checks on the SW parameters are performed
        if (Arrays.equals(refSeq, altSeq)) {
            final Cigar matching = new Cigar();
            matching.add(new CigarElement(refSeq.length, CigarOperator.MATCH_OR_MISMATCH));
            return matching;
        }

        final Cigar nonStandard;

        final String paddedRef = SW_PAD + new String(refSeq) + SW_PAD;
        final String paddedPath = SW_PAD + new String(altSeq) + SW_PAD;
        final SmithWatermanAlignment alignment = aligner.align(paddedRef.getBytes(), paddedPath.getBytes(), haplotypeToReferenceSWParameters, strategy);

        if ( isSWFailure(alignment) ) {
            return null;
        }

        // cut off the padding bases
        final int baseStart = SW_PAD.length();
        final int baseEnd = paddedPath.length() - SW_PAD.length() - 1; // -1 because it's inclusive
        final CigarBuilder.Result trimmedCigarAndDeletionsRemoved = AlignmentUtils.trimCigarByBases(alignment.getCigar(), baseStart, baseEnd);

        nonStandard = trimmedCigarAndDeletionsRemoved.getCigar();

        // leading deletion removed by cigar trimming shift the alignment start to the right
        final int trimmedLeadingDeletions = trimmedCigarAndDeletionsRemoved.getLeadingDeletionBasesRemoved();

        // trailing deletions should be kept in order to left-align
        final int trimmedTrailingDeletions = trimmedCigarAndDeletionsRemoved.getTrailingDeletionBasesRemoved();

        if ( trimmedTrailingDeletions > 0  ) {
            nonStandard.add(new CigarElement(trimmedTrailingDeletions, CigarOperator.D));
        }

        final CigarBuilder.Result leftAlignmentResult = AlignmentUtils.leftAlignIndels(nonStandard, refSeq, altSeq, trimmedLeadingDeletions);

        // we must account for possible leading deletions removed when trimming the padding and when left-aligning
        // trailing deletions removed when trimming have already been restored for left-alignment, but left-alingment may have removed them again.
        final int totalLeadingDeletionsRemoved = trimmedLeadingDeletions + leftAlignmentResult.getLeadingDeletionBasesRemoved();
        final int totalTrailingDeletionsRemoved = leftAlignmentResult.getTrailingDeletionBasesRemoved();

        if (totalLeadingDeletionsRemoved == 0 && totalTrailingDeletionsRemoved == 0) {
            return leftAlignmentResult.getCigar();
        } else {
            final List<CigarElement> resultElements = new ArrayList<>();
            if (totalLeadingDeletionsRemoved > 0) {
                resultElements.add(new CigarElement(totalLeadingDeletionsRemoved, CigarOperator.D));
            }
            resultElements.addAll(leftAlignmentResult.getCigar().getCigarElements());
            if (totalTrailingDeletionsRemoved > 0) {
                resultElements.add(new CigarElement(totalTrailingDeletionsRemoved, CigarOperator.D));
            }
            return new Cigar(resultElements);
        }
    }

    /**
     * Make sure that the SW didn't fail in some terrible way, and throw exception if it did
     */
    private static boolean isSWFailure(final SmithWatermanAlignment alignment) {
        // check that the alignment starts at the first base, which it should given the padding
        if ( alignment.getAlignmentOffset() > 0 ) {
            return true;
        }

        // check that we aren't getting any S operators (which would be very bad downstream)
        for ( final CigarElement ce : alignment.getCigar().getCigarElements() ) {
            if ( ce.getOperator() == CigarOperator.S )
                return true;
            // soft clips at the end of the alignment are really insertions
        }

        return false;
    }

    /**
     * Returns the length of the original read considering all clippings based on this cigar.
     * <p>
     *     The result of applying this method on a empty cigar is zero.
     * </p>
     * @param cigar the input cigar.
     * @throws IllegalArgumentException if the input {@code cigar} is {@code null}.
     * @return 0 or greater.
     */
    public static int countUnclippedReadBases(final Cigar cigar) {
        Utils.nonNull(cigar, "the input cigar cannot be null");
        return cigar.getCigarElements().stream()
                .filter(ce -> {
                    final CigarOperator operator = ce.getOperator();
                    return operator.isClipping() || operator.consumesReadBases();
                })
                .mapToInt(CigarElement::getLength)
                .sum();
    }

    /**
     * Count the number of soft- or hard- clipped bases from either the left or right end of a cigar
     */
    public static int countClippedBases(final Cigar cigar, final Tail tail, final CigarOperator typeOfClip) {
        Utils.validateArg(typeOfClip.isClipping(), "typeOfClip must be a clipping operator");

        final int size = cigar.numCigarElements();
        if (size < 2) {
            Utils.validateArg(size == 1 && !cigar.getFirstCigarElement().getOperator().isClipping(), "cigar is empty or completely clipped.");
            return 0;
        }

        int result = 0;

        for (int n = 0; n < size; n++) {
            final int index = (tail == Tail.LEFT ? n : size - n - 1);
            final CigarElement element = cigar.getCigarElement(index);
            if (!element.getOperator().isClipping()) {
                return result;
            } else if (element.getOperator() == typeOfClip) {
                result += element.getLength();
            }
        }

        throw new IllegalArgumentException("Input cigar " + cigar + " is completely clipped.");
    }

    /**
     * Count the number clipped bases (both soft and hard) from either the left or right end of a cigar
     */
    public static int countClippedBases(final Cigar cigar, final Tail tail) {
        return countClippedBases(cigar, tail, CigarOperator.SOFT_CLIP) + countClippedBases(cigar, tail, CigarOperator.HARD_CLIP);
    }

    /**
     * Count the number of soft- and hard-clipped bases over both ends of a cigar
     */
    public static int countClippedBases(final Cigar cigar, final CigarOperator clippingType) {
        return countClippedBases(cigar, Tail.LEFT, clippingType) + countClippedBases(cigar, Tail.RIGHT, clippingType);
    }

    /**
     * Count the number of clipped bases (both soft and hard) over both ends of a cigar
     */
    public static int countClippedBases(final Cigar cigar) {
        return countClippedBases(cigar, Tail.LEFT) + countClippedBases(cigar, Tail.RIGHT);
    }

    public static int countAlignedBases(final Cigar cigar ) {
        return Utils.nonNull(cigar).getCigarElements().stream()
                .filter(cigarElement -> cigarElement.getOperator().isAlignment())
                .mapToInt(CigarElement::getLength)
                .sum();
    }

    /**
     * replace soft clips (S) with match (M) operators, normalizing the result by all the transformations of the {@link CigarBuilder} class:
     * merging consecutive identical operators and removing zero-length elements.  For example 10S10M -> 20M and 10S10M10I10I -> 20M20I.
     */
    public static Cigar revertSoftClips(final Cigar originalCigar) {
        final CigarBuilder builder = new CigarBuilder();
        for (final CigarElement element : originalCigar.getCigarElements()) {
            if (element.getOperator() == CigarOperator.SOFT_CLIP) {
                builder.add(new CigarElement(element.getLength(), CigarOperator.MATCH_OR_MISMATCH));
            } else {
                builder.add(element);
            }
        }

        return builder.make();
    }

    /**
     * Given a cigar string, soft clip up to leftClipEnd and soft clip starting at rightClipBegin
     * @param start initial index to clip within read bases, inclusive
     * @param stop final index to clip within read bases exclusive
     * @param clippingOperator      type of clipping -- must be either hard clip or soft clip
     */
    public static Cigar clipCigar(final Cigar cigar, final int start, final int stop, CigarOperator clippingOperator) {
        Utils.validateArg(clippingOperator.isClipping(), "Not a clipping operator");
        final boolean clipLeft = start == 0;

        final CigarBuilder newCigar = new CigarBuilder();

        int elementStart = 0;
        for (final CigarElement element : cigar.getCigarElements()) {
            final CigarOperator operator = element.getOperator();
            // copy hard clips
            if (operator == CigarOperator.HARD_CLIP) {
                newCigar.add(new CigarElement(element.getLength(), element.getOperator()));
                continue;
            }
            final int elementEnd = elementStart + (operator.consumesReadBases() ? element.getLength() : 0);

            // element precedes start or follows end of clip, copy it to new cigar
            if (elementEnd <= start || elementStart >= stop) {
                // edge case: deletions at edge of clipping are meaningless and we skip them
                if (operator.consumesReadBases() || (elementStart != start && elementStart != stop)) {
                    newCigar.add(new CigarElement(element.getLength(), operator));
                }
            } else {    // otherwise, some or all of the element is soft-clipped
                final int unclippedLength = clipLeft ? elementEnd - stop : start - elementStart;
                final int clippedLength = element.getLength() - unclippedLength;

                if (unclippedLength <= 0) { // totally clipped
                    if (operator.consumesReadBases()) {
                        newCigar.add(new CigarElement(element.getLength(), clippingOperator));
                    }
                } else if (clipLeft) {
                    newCigar.add(new CigarElement(clippedLength, clippingOperator));
                    newCigar.add(new CigarElement(unclippedLength, operator));
                } else {
                    newCigar.add(new CigarElement(unclippedLength, operator));
                    newCigar.add(new CigarElement(clippedLength, clippingOperator));
                }
            }
            elementStart = elementEnd;
        }

        return newCigar.make();
    }

    /**
     * How many bases to the right does a read's alignment start shift given its cigar and the number of left soft clips
     */
    public static int alignmentStartShift(final Cigar cigar, final int numClipped) {
        int refBasesClipped = 0;

        int elementStart = 0;   // this and elementEnd are indices in the read's bases
        for (final CigarElement element : cigar.getCigarElements()) {
            final CigarOperator operator = element.getOperator();
            // hard clips consume neither read bases nor reference bases and are irrelevant
            if (operator == CigarOperator.HARD_CLIP) {
                continue;
            }
            final int elementEnd = elementStart + (operator.consumesReadBases() ? element.getLength() : 0);

            if (elementEnd <= numClipped) {  // totally within clipped span -- this includes deletions immediately following clipping
                refBasesClipped += operator.consumesReferenceBases() ? element.getLength() : 0;
            } else if (elementStart < numClipped) { // clip in middle of element, which means the element necessarily consumes read bases
                final int clippedLength = numClipped - elementStart;
                refBasesClipped += operator.consumesReferenceBases() ? clippedLength : 0;
                break;
            }
            elementStart = elementEnd;
        }
        return refBasesClipped;
    }

    /**
     * Computes the corresponding distance needs to be walked on the read, given the Cigar and distance walked on the reference.
     * @param cigar   cigar along the 5-3 direction of read (when read is mapped to reverse strand, bwa mem output cigar should be inverted)
     * @param start      start position (1-based) on the read (note it should not count the hard clipped bases, as usual)
     * @param refDist                 distance to walk on the reference
     * @param backward              whether to walk backwards along the read or not
     * @return                          corresponding walk distance on read (always positive)
     * @throws IllegalArgumentException if input cigar contains padding operation or 'N', or
     *                                  either of {@code start} or distance is non-positive, or
     *                                  {@code start} is larger than read length, or
     *                                  requested reference walk distance is longer than the total read bases in cigar, or
     *                                  computed read walk distance would "walk off" the read
     */
    public static int computeAssociatedDistOnRead(final Cigar cigar, final int start, final int refDist, final boolean backward) {

        Utils.validateArg(refDist > 0 && start > 0, () -> "start " + start + " or distance " + refDist + " is non-positive.");

        final List<CigarElement> elements = backward ? Lists.reverse(cigar.getCigarElements()) : cigar.getCigarElements();

        final int readLength = elements.stream().mapToInt(ce -> ce.getOperator().consumesReadBases() ? ce.getLength() : 0).sum();
        final int readBasesToSkip = backward ? readLength - start : start - 1;

        int readBasesConsumed = 0;
        int refBasesConsumed = 0;

        for (final CigarElement element : elements){
            final int readBasesConsumedBeforeElement = readBasesConsumed;

            readBasesConsumed += element.getOperator().consumesReadBases() ? element.getLength() : 0;
            // skip cigar elements that end before the read start or start after the reference end
            if (readBasesConsumed <= readBasesToSkip) {
                continue;
            }

            refBasesConsumed += element.getOperator().consumesReferenceBases() ? element.getLength() - Math.max(readBasesToSkip - readBasesConsumedBeforeElement, 0) : 0;
            if (refBasesConsumed >= refDist) {
                final int excessRefBasesInElement = Math.max(refBasesConsumed - refDist, 0);
                return readBasesConsumed - readBasesToSkip - (element.getOperator().consumesReadBases() ? excessRefBasesInElement : 0);
            }
        }

        throw new IllegalArgumentException("Cigar " + cigar + "does not contain at least " + refDist + " reference bases past red start " + start + ".");
    }

    /**
     * Convert the 'I' CigarElement, if it is at either end (terminal) of the input cigar, to a corresponding 'S' operator.
     * Note that we allow CIGAR of the format '10H10S10I10M', but disallows the format if after the conversion the cigar turns into a giant clip,
     * e.g. '10H10S10I10S10H' is not allowed (if allowed, it becomes a giant clip of '10H30S10H' which is non-sense).
     *
     * @return a pair of number of clipped (hard and soft, including the ones from the converted terminal 'I') bases at the front and back of the
     *         input {@code cigarAlongInput5to3Direction}.
     *
     * @throws IllegalArgumentException when the checks as described above fail.
     */
    public static Cigar convertTerminalInsertionToSoftClip(final Cigar cigar) {

        if (cigar.numCigarElements() < 2 ) {
            return cigar;
        }

        final CigarBuilder builder = new CigarBuilder();
        for (int n = 0; n < cigar.numCigarElements(); n++) {
            final CigarElement element = cigar.getCigarElement(n);
            if (element.getOperator() != CigarOperator.INSERTION) { // not an insertion
                builder.add(element);
            } else if (n == 0 || n == cigar.numCigarElements() - 1) {   // terminal insertion with no clipping -- convert to soft clip
                builder.add(new CigarElement(element.getLength(), CigarOperator.SOFT_CLIP));
            } else if (cigar.getCigarElement(n-1).getOperator().isClipping() || cigar.getCigarElement(n+1).getOperator().isClipping()) {    // insertion preceding or following clip
                builder.add(new CigarElement(element.getLength(), CigarOperator.SOFT_CLIP));
            } else {    // interior insertion
                builder.add(element);
            }
        }

        return builder.make();
    }
}
