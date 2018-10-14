package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;

import java.util.*;

public final class CigarUtils {

    // used in the bubble state machine to apply Smith-Waterman to the bubble sequence
    // these values were chosen via optimization against the NA12878 knowledge base
    public static final SWParameters NEW_SW_PARAMETERS = new SWParameters(200, -150, -260, -11);

    // In Mutect2 and HaplotypeCaller reads are realigned to their *best* haplotypes, which is very different from a generic alignment.
    // The {@code NEW_SW_PARAMETERS} penalize a substitution error more than an indel up to a length of 9 bases!
    // Suppose, for example, that a read has a single substitution error, say C -> T, on its last base.  Those parameters
    // would prefer to extend a deletion until the next T on the reference is found in order to avoid the substitution, which is absurd.
    // Since these parameters are for aligning a read to the biological sequence we believe it comes from, the parameters
    // we choose should correspond to sequencer error.  They *do not* have anything to do with the prevalence of true variation!
    public static final SWParameters ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS = new SWParameters(10, -15, -30, -5);

    private static final String SW_PAD = "NNNNNNNNNN";

    private CigarUtils(){}

    /**
     * Combines equal adjacent elements of a Cigar object
     *
     * @param rawCigar the cigar object
     * @return a combined cigar object
     */
    public static Cigar combineAdjacentCigarElements(final Cigar rawCigar) {
        Utils.nonNull(rawCigar);
        final Cigar combinedCigar = new Cigar();
        CigarElement lastElement = null;
        int lastElementLength = 0;
        for (final CigarElement cigarElement : rawCigar.getCigarElements()) {
            if (lastElement != null &&
                    ((lastElement.getOperator() == cigarElement.getOperator()) ||
                            (lastElement.getOperator() == CigarOperator.I && cigarElement.getOperator() == CigarOperator.D) ||
                            (lastElement.getOperator() == CigarOperator.D && cigarElement.getOperator() == CigarOperator.I)))
                lastElementLength += cigarElement.getLength();
            else
            {
                if (lastElement != null)
                    combinedCigar.add(new CigarElement(lastElementLength, lastElement.getOperator()));

                lastElement = cigarElement;
                lastElementLength = cigarElement.getLength();
            }
        }
        if (lastElement != null) {
            combinedCigar.add(new CigarElement(lastElementLength, lastElement.getOperator()));
        }

        return combinedCigar;
    }

    /**
     * Checks whether the cigar has any element that is not H or S
     * @return true the cigar has elements other than S or H, false otherwise.
     */
    public static boolean hasNonClippedBases(final Cigar cigar) {
        return Utils.nonNull(cigar).getCigarElements().stream()
                .anyMatch(el -> el.getOperator() != CigarOperator.SOFT_CLIP && el.getOperator() != CigarOperator.HARD_CLIP);
    }

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
    public static int countRefBasesBasedOnUnclippedAlignment(final GATKRead read, final int cigarStartIndex, final int cigarEndIndex){
        if (read == null){
            throw new IllegalArgumentException("null read");
        }
        final List<CigarElement> elems = read.getCigarElements();
        if (cigarStartIndex < 0 || cigarEndIndex > elems.size() || cigarStartIndex > cigarEndIndex){
            throw new IllegalArgumentException("invalid index:" + 0 + " -" + elems.size());
        }
        int result = 0;
        for(int i = cigarStartIndex; i < cigarEndIndex; i++){
            final CigarElement cigarElement = elems.get(i);
            final CigarOperator operator = cigarElement.getOperator();
            if (operator.consumesReferenceBases() || operator.isClipping()) {
                result += cigarElement.getLength();
            }
        }
        return result;
    }

    public static int countRefBasesIncludingSoftClips(final GATKRead read, final int cigarStartIndex, final int cigarEndIndex){
        Utils.nonNull(read, "null read");
        final List<CigarElement> elems = read.getCigarElements();
        if (cigarStartIndex < 0 || cigarEndIndex > elems.size() || cigarStartIndex > cigarEndIndex){
            throw new IllegalArgumentException("invalid index:" + 0 + " -" + elems.size());
        }
        int result = 0;
        for(int i = cigarStartIndex; i < cigarEndIndex; i++){
            final CigarElement cigarElement = elems.get(i);
            final CigarOperator operator = cigarElement.getOperator();
            if (operator.consumesReferenceBases() || operator == CigarOperator.S) {
                result += cigarElement.getLength();
            }
        }
        return result;
    }

    /**
     * Removes all clipping operators from the cigar.
     */
    public static Cigar trimReadToUnclippedBases(final Cigar cigar) {
        Utils.nonNull(cigar, "cigar is null");
        final List<CigarElement> elements = new ArrayList<>(cigar.numCigarElements());
        for ( final CigarElement ce : cigar.getCigarElements() ) {
            if ( !isClipOperator(ce.getOperator()) ) {
                elements.add(ce);
            }
        }
        return new Cigar(elements);
    }

    private static boolean isClipOperator(final CigarOperator op) {
        return op == CigarOperator.S || op == CigarOperator.H || op == CigarOperator.P;
    }

    private static boolean isClipOperator(final CigarElement el) {
        return isClipOperator(el.getOperator());
    }

    /**
     * Given a cigar1 and a read with cigar2,
     * this method creates cigar3 such that it has flanking clip operators from cigar2
     * and it has all operators from cigar1 in the middle.
     *
     * In other words if:
     * cigar2 = leftClip2 + noclips2 + rightClip2
     *
     * then
     * cigar3 = leftClip2 + cigar1 + rightClip2
     */
    public static Cigar reclipCigar(final Cigar cigar, final GATKRead read) {
        Utils.nonNull(cigar, "cigar");
        Utils.nonNull(read, "read");

        final List<CigarElement> elements = new ArrayList<>();
        int i = 0;
        final Cigar readCigar = read.getCigar();
        final int n = readCigar.numCigarElements();
        final List<CigarElement> readEls = readCigar.getCigarElements();

        //copy head clips
        while ( i < n && isClipOperator(readEls.get(i)) ) {
            elements.add(readEls.get(i));
            i++;
        }

        elements.addAll(cigar.getCigarElements());

        //skip over non-clips
        i++;
        while ( i < n && !isClipOperator(readEls.get(i)) ) {
            i++;
        }

        //copy tail clips
        while ( i < n && isClipOperator(readEls.get(i)) ) {
            elements.add(readEls.get(i));
            i++;
        }

        return new Cigar(elements);
    }

    /**
     * Returns whether the cigar has any N operators.
     */
    public static boolean containsNOperator(final Cigar cigar) {
        Utils.nonNull(cigar);
        //Note: reach the elements directly rather that calling getCigarElements because
        // we want to avoid allocating a new unmodifiable list view (comes up in profiling of HaplotypeCaller)
        for (int i = 0, n = cigar.numCigarElements(); i < n; i++) {
            if (cigar.getCigarElement(i).getOperator() == CigarOperator.N){
                return true;
            }
        }
        return false;
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
        if (hasConsecutiveIndels(elems)){
            return false;
        }
        if (startsWithDeletionIgnoringClips(elems)){
            return false;
        }
        //revert the list and check deletions at the end
        final List<CigarElement> elemsRev = new ArrayList<>(elems);
        Collections.reverse(elemsRev);
        return !startsWithDeletionIgnoringClips(elemsRev);
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
    private static boolean startsWithDeletionIgnoringClips(final List<CigarElement> elems) {
        final Iterator<CigarElement> iter = elems.iterator();
        boolean isClip = true;
        CigarOperator op = null;
        while(iter.hasNext() && isClip) { //consume clips at the beginning
            final CigarElement elem = iter.next();
            op = elem.getOperator();
            isClip = (op == CigarOperator.HARD_CLIP || op == CigarOperator.SOFT_CLIP);
        }
        //once all clips are consumed, is it a deletion or not?
        return op == CigarOperator.DELETION;
    }

    /**
     * Calculate the cigar elements for this path against the reference sequence
     *
     * @param aligner
     * @param refSeq the reference sequence that all of the bases in this path should align to
     * @return a Cigar mapping this path to refSeq, or null if no reasonable alignment could be found
     */
    public static Cigar calculateCigar(final byte[] refSeq, final byte[] altSeq, final SmithWatermanAligner aligner) {
        Utils.nonNull(refSeq, "refSeq");
        Utils.nonNull(altSeq, "altSeq");
        if ( altSeq.length == 0 ) {
            // horrible edge case from the unit tests, where this path has no bases
            return new Cigar(Collections.singletonList(new CigarElement(refSeq.length, CigarOperator.D)));
        }

        //Note: this is a performance optimization.
        // If two strings are equal (a O(n) check) then it's trivial to get CIGAR for them.
        if (Arrays.equals(refSeq, altSeq)){
            final Cigar matching = new Cigar();
            matching.add(new CigarElement(refSeq.length, CigarOperator.MATCH_OR_MISMATCH));
            return matching;
        }

        final Cigar nonStandard;

        final String paddedRef = SW_PAD + new String(refSeq) + SW_PAD;
        final String paddedPath = SW_PAD + new String(altSeq) + SW_PAD;
        final SmithWatermanAlignment alignment = aligner.align(paddedRef.getBytes(), paddedPath.getBytes(), NEW_SW_PARAMETERS, SWOverhangStrategy.SOFTCLIP);

        if ( isSWFailure(alignment) ) {
            return null;
        }


        // cut off the padding bases
        final int baseStart = SW_PAD.length();
        final int baseEnd = paddedPath.length() - SW_PAD.length() - 1; // -1 because it's inclusive
        nonStandard = AlignmentUtils.trimCigarByBases(alignment.getCigar(), baseStart, baseEnd);

        if ( nonStandard.getReferenceLength() != refSeq.length ) {
            nonStandard.add(new CigarElement(refSeq.length - nonStandard.getReferenceLength(), CigarOperator.D));
        }

        // finally, return the cigar with all indels left aligned
        return leftAlignCigarSequentially(nonStandard, refSeq, altSeq, 0, 0);
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
     * Left align the given cigar sequentially. This is needed because AlignmentUtils doesn't accept cigars with more than one indel in them.
     * This is a target of future work to incorporate and generalize into AlignmentUtils for use by others.
     * @param cigar     the cigar to left align
     * @param refSeq    the reference byte array
     * @param readSeq   the read byte array
     * @param refIndex  0-based alignment start position on ref
     * @param readIndex 0-based alignment start position on read
     * @return          the left-aligned cigar
     */
    public static Cigar leftAlignCigarSequentially(final Cigar cigar, final byte[] refSeq, final byte[] readSeq, int refIndex, int readIndex) {
        Utils.nonNull(cigar, "cigar null");
        Utils.nonNull(refSeq, "refSeq null");
        Utils.nonNull(readSeq, "readSeq null");

        final Cigar cigarToReturn = new Cigar();
        Cigar cigarToAlign = new Cigar();
        for (int i = 0; i < cigar.numCigarElements(); i++) {
            final CigarElement ce = cigar.getCigarElement(i);
            if (ce.getOperator() == CigarOperator.D || ce.getOperator() == CigarOperator.I) {
                cigarToAlign.add(ce);
                final Cigar leftAligned = AlignmentUtils.leftAlignSingleIndel(cigarToAlign, refSeq, readSeq, refIndex, readIndex, false);
                for ( final CigarElement toAdd : leftAligned.getCigarElements() ) { cigarToReturn.add(toAdd); }
                refIndex += cigarToAlign.getReferenceLength();
                readIndex += cigarToAlign.getReadLength();
                cigarToAlign = new Cigar();
            } else {
                cigarToAlign.add(ce);
            }
        }
        if( !cigarToAlign.isEmpty() ) {
            for( final CigarElement toAdd : cigarToAlign.getCigarElements() ) {
                cigarToReturn.add(toAdd);
            }
        }

        final Cigar result = AlignmentUtils.consolidateCigar(cigarToReturn);
        Utils.validate(result.getReferenceLength() == cigar.getReferenceLength(),
                () -> "leftAlignCigarSequentially failed to produce a valid CIGAR.  Reference lengths differ.  Initial cigar " + cigar + " left aligned into " + result);
        return result;
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
     * Total number of bases clipped on the left/head side of the cigar.
     *
     * @param cigar the input cigar.
     * @throws IllegalArgumentException if {@code cigar} is {@code null}.
     * @return 0 or greater.
     */
    public static int countLeftClippedBases(final Cigar cigar) {
        Utils.nonNull(cigar, "the input cigar cannot not be null");
        if (cigar.numCigarElements() < 2) {
            return 0;
        } else {
            int result = 0;
            for (final CigarElement e : cigar) {
                if (!e.getOperator().isClipping()) {
                    return result;
                }
                result += e.getLength();
            }
            throw new IllegalArgumentException("the input cigar only contains clips!");
        }
    }

    /**
     * Returns the number of based hard-clipped to the left/head of the cigar.
     *
     * @param cigar the input cigar.
     * @throws IllegalArgumentException if {@code cigar} is {@code null}.
     * @return 0 or greater.
     */
    public static int countLeftHardClippedBases(final Cigar cigar) {
        Utils.nonNull(cigar, "the input cigar cannot not be null");
        if (cigar.numCigarElements() < 2) {
            return 0;
        } else if (cigar.getCigarElement(0).getOperator() != CigarOperator.H) {
            return 0;
        } else {
            return cigar.getCigarElement(0).getLength();
        }
    }

    /**
     * Returns the number of based hard-clipped to the right/tail of the cigar.
     *
     * @param cigar the input cigar.
     * @throws IllegalArgumentException if {@code cigar} is {@code null}.
     * @return 0 or greater.
     */
    public static int countRightHardClippedBases(final Cigar cigar) {
        Utils.nonNull(cigar, "the input cigar cannot not be null");
        if (cigar.numCigarElements() < 2) {
            return 0;
        } else {
            final List<CigarElement> elements = cigar.getCigarElements();
            int lastElementIndex;
            if (elements.get(lastElementIndex = elements.size() - 1).getOperator() != CigarOperator.H) {
                return 0;
            } else {
                return elements.get(lastElementIndex).getLength();
            }
        }
    }

    /**
     * Total number of bases clipped (soft or hard) on the right/tail side of the cigar.
     *
     * @param cigar the input cigar.
     * @throws IllegalArgumentException if {@code cigar} is {@code null}
     * @return 0 or greater.
     */
    public static int countRightClippedBases(final Cigar cigar) {
        Utils.nonNull(cigar, "the input cigar cannot be null");
        final List<CigarElement> elements = cigar.getCigarElements();
        final int elementsCount = elements.size();
        if (elementsCount < 2) {  // a single clipped element (that is already an "invalid" CIGAR) would be considered
            return 0;        // a left-clip so it must have at least two elements before right clipping can be larger than 0.
        } else {
            int result = 0;
            int i;
            for (i = elementsCount - 1; i >= 0; --i) {
                final CigarElement ce = elements.get(i);
                if (!ce.getOperator().isClipping()) {
                    return result;
                } else {
                    result += ce.getLength();
                }
            }
            throw new IllegalArgumentException("the input cigar only have clipping operations");
        }
    }

    public static int countAlignedBases(final Cigar cigar ) {
        return Utils.nonNull(cigar).getCigarElements().stream()
                .filter(cigarElement -> cigarElement.getOperator().isAlignment())
                .mapToInt(CigarElement::getLength)
                .sum();
    }
}
