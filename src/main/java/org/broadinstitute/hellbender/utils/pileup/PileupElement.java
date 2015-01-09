/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.utils.pileup;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.sam.ReadUtils;

import java.util.*;

public class PileupElement implements Comparable<PileupElement> {
    private final static LinkedList<CigarElement> EMPTY_LINKED_LIST = new LinkedList<>();

    private final static EnumSet<CigarOperator> ON_GENOME_OPERATORS =
            EnumSet.of(CigarOperator.M, CigarOperator.EQ, CigarOperator.X, CigarOperator.D);

    public static final byte DELETION_BASE = BaseUtils.Base.D.base;
    public static final byte DELETION_QUAL = (byte) 16;
    public static final byte A_FOLLOWED_BY_INSERTION_BASE = (byte) 87;
    public static final byte C_FOLLOWED_BY_INSERTION_BASE = (byte) 88;
    public static final byte T_FOLLOWED_BY_INSERTION_BASE = (byte) 89;
    public static final byte G_FOLLOWED_BY_INSERTION_BASE = (byte) 90;

    protected final SAMRecord read;         // the read this base belongs to
    protected final int offset;                 // the offset in the bases array for this base

    private final CigarElement currentCigarElement;
    private final int currentCigarOffset;
    private final int offsetInCurrentCigar;

    /**
     * Create a new pileup element
     *
     * @param read a non-null read to pileup
     * @param baseOffset the offset into the read's base / qual vector aligned to this position on the genome. If the
     *                   current cigar element is a deletion, offset should be the offset of the last M/=/X position.
     * @param currentElement a non-null CigarElement that indicates the cigar element aligning the read to the genome
     * @param currentCigarOffset the offset of currentElement in read.getCigar().getElement(currentCigarOffset) == currentElement)
     * @param offsetInCurrentCigar how far into the currentElement are we in our alignment to the genome?
     */
    public PileupElement(final SAMRecord read, final int baseOffset,
                         final CigarElement currentElement, final int currentCigarOffset,
                         final int offsetInCurrentCigar) {
        assert currentElement != null;

        this.read = read;
        this.offset = baseOffset;
        this.currentCigarElement = currentElement;
        this.currentCigarOffset = currentCigarOffset;
        this.offsetInCurrentCigar = offsetInCurrentCigar;

        // for performance regions these are assertions
        assert this.read != null;
        assert this.offset >= 0 && this.offset < this.read.getReadLength();
        assert this.currentCigarOffset >= 0;
        assert this.currentCigarOffset < read.getCigarLength();
        assert this.offsetInCurrentCigar >= 0;
        assert this.offsetInCurrentCigar < currentElement.getLength();
    }

    /**
     * Create a new PileupElement that's a copy of toCopy
     * @param toCopy the element we want to copy
     */
    public PileupElement(final PileupElement toCopy) {
        this(toCopy.read, toCopy.offset, toCopy.currentCigarElement, toCopy.currentCigarOffset, toCopy.offsetInCurrentCigar);
    }

    /**
     * Is this element a deletion w.r.t. the reference genome?
     *
     * @return true if this is a deletion, false otherwise
     */
    public boolean isDeletion() {
        return currentCigarElement.getOperator() == CigarOperator.D;
    }

    /**
     * Is the current element immediately before a deletion, but itself not a deletion?
     *
     * Suppose we are aligning a read with cigar 3M2D1M.  This function is true
     * if we are in the last cigar position of the 3M, but not if we are in the 2D itself.
     *
     * @return true if the next alignment position is a deletion w.r.t. the reference genome
     */
    public boolean isBeforeDeletionStart() {
        return ! isDeletion() && atEndOfCurrentCigar() && hasOperator(getNextOnGenomeCigarElement(), CigarOperator.D);
    }

    /**
     * Is the current element immediately after a deletion, but itself not a deletion?
     *
     * Suppose we are aligning a read with cigar 1M2D3M.  This function is true
     * if we are in the first cigar position of the 3M, but not if we are in the 2D itself or
     * in any but the first position of the 3M.
     *
     * @return true if the previous alignment position is a deletion w.r.t. the reference genome
     */
    public boolean isAfterDeletionEnd() {
        return ! isDeletion() && atStartOfCurrentCigar() && hasOperator(getPreviousOnGenomeCigarElement(), CigarOperator.D);
    }

    /**
     * Get the read for this pileup element
     * @return a non-null GATKSAMRecord
     */
    public SAMRecord getRead() {
        return read;
    }

    /**
     * Get the offset of the this element into the read that aligns that read's base to this genomic position.
     *
     * If the current element is a deletion then offset is the offset of the last base containing offset.
     *
     * @return a valid offset into the read's bases
     */
    public int getOffset() {
        return offset;
    }

    /**
     * Get the base aligned to the genome at this location
     *
     * If the current element is a deletion returns DELETION_BASE
     *
     * @return a base encoded as a byte
     */
    public byte getBase() {
        return isDeletion() ? DELETION_BASE : read.getReadBases()[offset];
    }

    /**
     * Get the base quality score of the base at this aligned position on the genome
     * @return a phred-scaled quality score as a byte
     */
    public byte getQual() {
        return isDeletion() ? DELETION_QUAL : read.getBaseQualities()[offset];
    }

    /**
     * Get the Base Insertion quality at this pileup position
     * @return a phred-scaled quality score as a byte
     */
    public byte getBaseInsertionQual() {
        return isDeletion() ? DELETION_QUAL : ReadUtils.getBaseInsertionQualities(read)[offset];
    }

    /**
     * Get the Base Deletion quality at this pileup position
     * @return a phred-scaled quality score as a byte
     */
    public byte getBaseDeletionQual() {
        return isDeletion() ? DELETION_QUAL : ReadUtils.getBaseDeletionQualities(read)[offset];
    }

    /**
     * Get the length of an immediately following insertion or deletion event, or 0 if no such event exists
     *
     * Only returns a positive value when this pileup element is immediately before an indel.  Being
     * immediately before a deletion means that this pileup element isn't an deletion, and that the
     * next genomic alignment for this read is a deletion.  For the insertion case, this means
     * that an insertion cigar occurs immediately after this element, between this one and the
     * next genomic position.
     *
     * Note this function may be expensive, so multiple uses should be cached by the caller
     *
     * @return length of the event (number of inserted or deleted bases), or 0
     */
    public int getLengthOfImmediatelyFollowingIndel() {
        final CigarElement element = getNextIndelCigarElement();
        return element == null ? 0 : element.getLength();
    }

    /**
     * Helpful function to get the immediately following cigar element, for an insertion or deletion
     *
     * if this state precedes a deletion (i.e., next position on genome) or insertion (immediately between
     * this and the next position) returns the CigarElement corresponding to this event.  Otherwise returns
     * null.
     *
     * @return a CigarElement, or null if the next alignment state ins't an insertion or deletion.
     */
    private CigarElement getNextIndelCigarElement() {
        if ( isBeforeDeletionStart() ) {
            final CigarElement element = getNextOnGenomeCigarElement();
            if ( element == null || element.getOperator() != CigarOperator.D )
                throw new IllegalStateException("Immediately before deletion but the next cigar element isn't a deletion " + element);
            return element;
        } else if ( isBeforeInsertion() ) {
            final CigarElement element = getBetweenNextPosition().get(0);
            if ( element.getOperator() != CigarOperator.I )
                throw new IllegalStateException("Immediately before insertion but the next cigar element isn't an insertion " + element);
            return element;
        } else {
            return null;
        }
    }

    /**
     * Get the bases for an insertion that immediately follows this alignment state, or null if none exists
     *
     * @see #getLengthOfImmediatelyFollowingIndel() for details on the meaning of immediately.
     *
     * If the immediately following state isn't an insertion, returns null
     *
     * @return actual sequence of inserted bases, or a null if the event is a deletion or if there is no event in the associated read.
     */
    public String getBasesOfImmediatelyFollowingInsertion() {
        final CigarElement element = getNextIndelCigarElement();
        if ( element != null && element.getOperator() == CigarOperator.I ) {
            final int getFrom = offset + 1;
            final byte[] bases = Arrays.copyOfRange(read.getReadBases(), getFrom, getFrom + element.getLength());
            return new String(bases);
        } else
            return null;
    }

    /**
     * Get the mapping quality of the read of this element
     * @return the mapping quality of the underlying SAM record
     */
    public int getMappingQual() {
        return read.getMappingQuality();
    }

    public String toString() {
        return String.format("%s @ %d = %c Q%d", getRead().getReadName(), getOffset(), (char) getBase(), getQual());
    }

    @Override
    public int compareTo(final PileupElement pileupElement) {
        if (offset < pileupElement.offset)
            return -1;
        else if (offset > pileupElement.offset)
            return 1;
        else if (read.getAlignmentStart() < pileupElement.read.getAlignmentStart())
            return -1;
        else if (read.getAlignmentStart() > pileupElement.read.getAlignmentStart())
            return 1;
        else
            return 0;
    }

    // --------------------------------------------------------------------------
    //
    // Reduced read accessors
    //
    // --------------------------------------------------------------------------

    /**
     * Get the cigar element aligning this element to the genome
     * @return a non-null CigarElement
     */
    public CigarElement getCurrentCigarElement() {
        return currentCigarElement;
    }

    /**
     * Get the offset of this cigar element in the Cigar of the current read (0-based)
     *
     * Suppose the cigar is 1M2D3I4D.  If we are in the 1M state this function returns
     * 0.  If we are in 2D, the result is 1.  If we are in the 4D, the result is 3.
     *
     * @return an offset into the read.getCigar() that brings us to the current cigar element
     */
    public int getCurrentCigarOffset() {
        return currentCigarOffset;
    }

    /**
     * Get the offset into the *current* cigar element for this alignment position
     *
     * We can be anywhere from offset 0 (first position) to length - 1 of the current
     * cigar element aligning us to this genomic position.
     *
     * @return a valid offset into the current cigar element
     */
    public int getOffsetInCurrentCigar() {
        return offsetInCurrentCigar;
    }

    /**
     * Get the cigar elements that occur before the current position but after the previous position on the genome
     *
     * For example, if we are in the 3M state of 1M2I3M state then 2I occurs before this position.
     *
     * Note that this function does not care where we are in the current cigar element.  In the previous
     * example this list of elements contains the 2I state regardless of where you are in the 3M.
     *
     * Note this returns the list of all elements that occur between this and the prev site, so for
     * example we might have 5S10I2M and this function would return [5S, 10I].
     *
     * @return a non-null list of CigarElements
     */
    public LinkedList<CigarElement> getBetweenPrevPosition() {
        return atStartOfCurrentCigar() ? getBetween(Direction.PREV) : EMPTY_LINKED_LIST;
    }

    /**
     * Get the cigar elements that occur after the current position but before the next position on the genome
     *
     * @see #getBetweenPrevPosition() for more details
     *
     * @return a non-null list of CigarElements
     */
    public List<CigarElement> getBetweenNextPosition() {
        return atEndOfCurrentCigar() ? getBetween(Direction.NEXT) : EMPTY_LINKED_LIST;
    }

    /** for some helper functions */
    private enum Direction { PREV, NEXT }

    /**
     * Helper function to get cigar elements between this and either the prev or next genomic position
     *
     * @param direction PREVIOUS if we want before, NEXT if we want after
     * @return a non-null list of cigar elements between this and the neighboring position in direction
     */
    private LinkedList<CigarElement> getBetween(final Direction direction) {
        final int increment = direction == Direction.NEXT ? 1 : -1;
        LinkedList<CigarElement> elements = null;
        final int nCigarElements = read.getCigarLength();
        for ( int i = currentCigarOffset + increment; i >= 0 && i < nCigarElements; i += increment) {
            final CigarElement elt = read.getCigar().getCigarElement(i);
            if ( ON_GENOME_OPERATORS.contains(elt.getOperator()) )
                break;
            else {
                // optimization: don't allocate list if not necessary
                if ( elements == null )
                    elements = new LinkedList<CigarElement>();

                if ( increment > 0 )
                    // to keep the list in the right order, if we are incrementing positively add to the end
                    elements.add(elt);
                else
                    // counting down => add to front
                    elements.addFirst(elt);
            }
        }

        // optimization: elements is null because nothing got added, just return the empty list
        return elements == null ? EMPTY_LINKED_LIST : elements;
    }

    /**
     * Get the cigar element of the previous genomic aligned position
     *
     * For example, we might have 1M2I3M, and be sitting at the someone in the 3M.  This
     * function would return 1M, as the 2I isn't on the genome.  Note this function skips
     * all of the positions that would occur in the current element.  So the result
     * is always 1M regardless of whether we're in the first, second, or third position of the 3M
     * cigar.
     *
     * @return a CigarElement, or null (indicating that no previous element exists)
     */
    public CigarElement getPreviousOnGenomeCigarElement() {
        return getNeighboringOnGenomeCigarElement(Direction.PREV);
    }

    /**
     * Get the cigar element of the next genomic aligned position
     *
     * @see #getPreviousOnGenomeCigarElement() for more details
     *
     * @return a CigarElement, or null (indicating that no next element exists)
     */
    public CigarElement getNextOnGenomeCigarElement() {
        return getNeighboringOnGenomeCigarElement(Direction.NEXT);
    }

    /**
     * Helper function to get the cigar element of the next or previous genomic position
     * @param direction the direction to look in
     * @return a CigarElement, or null if no such element exists
     */
    private CigarElement getNeighboringOnGenomeCigarElement(final Direction direction) {
        final int increment = direction == Direction.NEXT ? 1 : -1;
        final int nCigarElements = read.getCigarLength();

        for ( int i = currentCigarOffset + increment; i >= 0 && i < nCigarElements; i += increment) {
            final CigarElement elt = read.getCigar().getCigarElement(i);
            if ( ON_GENOME_OPERATORS.contains(elt.getOperator()) )
                return elt;
        }

        // getting here means that you didn't find anything
        return null;
    }

    /**
     * Does the cigar element (which may be null) have operation toMatch?
     *
     * @param maybeCigarElement a CigarElement that might be null
     * @param toMatch a CigarOperator we want to match against the one in maybeCigarElement
     * @return true if maybeCigarElement isn't null and has operator toMatch
     */
    private boolean hasOperator(final CigarElement maybeCigarElement, final CigarOperator toMatch) {
        return maybeCigarElement != null && maybeCigarElement.getOperator() == toMatch;
    }

    /**
     * Does an insertion occur immediately before the current position on the genome?
     *
     * @return true if yes, false if no
     */
    public boolean isAfterInsertion() { return isAfter(getBetweenPrevPosition(), CigarOperator.I); }

    /**
     * Does an insertion occur immediately after the current position on the genome?
     *
     * @return true if yes, false if no
     */
    public boolean isBeforeInsertion() { return isBefore(getBetweenNextPosition(), CigarOperator.I); }

    /**
     * Does a soft-clipping event occur immediately before the current position on the genome?
     *
     * @return true if yes, false if no
     */
    public boolean isAfterSoftClip() { return isAfter(getBetweenPrevPosition(), CigarOperator.S); }

    /**
     * Does a soft-clipping event occur immediately after the current position on the genome?
     *
     * @return true if yes, false if no
     */
    public boolean isBeforeSoftClip() { return isBefore(getBetweenNextPosition(), CigarOperator.S); }

    /**
     * Does a soft-clipping event occur immediately before or after the current position on the genome?
     *
     * @return true if yes, false if no
     */
    public boolean isNextToSoftClip() { return isAfterSoftClip() || isBeforeSoftClip(); }

    /**
     * Is the current position at the end of the current cigar?
     *
     * For example, if we are in element 3M, this function returns true if we are at offsetInCurrentCigar
     * of 2, but not 0 or 1.
     *
     * @return true if we're at the end of the current cigar
     */
    public boolean atEndOfCurrentCigar() {
        return offsetInCurrentCigar == currentCigarElement.getLength() - 1;
    }

    /**
     * Is the current position at the start of the current cigar?
     *
     * For example, if we are in element 3M, this function returns true if we are at offsetInCurrentCigar
     * of 0, but not 1 or 2.
     *
     * @return true if we're at the start of the current cigar
     */
    public boolean atStartOfCurrentCigar() {
        return offsetInCurrentCigar == 0;
    }

    /**
     * Is op the last element in the list of elements?
     *
     * @param elements the elements to examine
     * @param op the op we want the last element's op to equal
     * @return true if op == last(elements).op
     */
    private boolean isAfter(final Deque<CigarElement> elements, final CigarOperator op) {
        return ! elements.isEmpty() && elements.peekLast().getOperator() == op;
    }

    /**
     * Is op the first element in the list of elements?
     *
     * @param elements the elements to examine
     * @param op the op we want the last element's op to equal
     * @return true if op == first(elements).op
     */
    private boolean isBefore(final List<CigarElement> elements, final CigarOperator op) {
        return ! elements.isEmpty() && elements.get(0).getOperator() == op;
    }
}
