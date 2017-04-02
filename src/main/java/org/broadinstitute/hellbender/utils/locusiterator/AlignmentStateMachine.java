package org.broadinstitute.hellbender.utils.locusiterator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Steps a single read along its alignment to the genome
 *
 * The logical model for generating extended events is as follows: the "record state"
 * implements the traversal along the reference; thus stepForwardOnGenome() returns
 * on every and only on actual reference bases. This can be a (mis)match or a deletion
 * (in the latter case, we still return on every individual reference base the deletion spans).
 *
 */
public final class AlignmentStateMachine {
    /**
     * Our read
     */
    private final GATKRead read;
    private final Cigar cigar;
    private final int nCigarElements;
    private int currentCigarElementOffset = -1;

    /**
     * how far are we offset from the start of the read bases?
     */
    private int readOffset;

    /**
     * how far are we offset from the alignment start on the genome?
     */
    private int genomeOffset;

    /**
     * Our cigar element
     */
    private CigarElement currentElement;

    /**
     * how far are we into our cigarElement?
     */
    private int offsetIntoCurrentCigarElement;

    public AlignmentStateMachine(final GATKRead read) {
        this.read = read;
        this.cigar = read.getCigar();
        this.nCigarElements = cigar.numCigarElements();
        initializeAsLeftEdge();
    }

    /**
     * Initialize the state variables to put this machine one bp before the
     * start of the alignment, so that a call to stepForwardOnGenome() will advance
     * us to the first proper location
     */
    private void initializeAsLeftEdge() {
        readOffset = offsetIntoCurrentCigarElement = genomeOffset = -1;
        currentElement = null;
    }

    /**
     * Get the read we are aligning to the genome
     * @return a non-null Read
     */
    public GATKRead getRead() {
        return read;
    }

    /**
     * Get the contig of the underlying read
     *
     * @return the contig of the read
     */
    public String getContig() {
        return getRead().getContig();
    }

    /**
     * Is this the left edge state?  I.e., one that is before or after the current read?
     * @return true if this state is an edge state, false otherwise
     */
    public boolean isLeftEdge() {
        return readOffset == -1;
    }

    /**
     * Are we on the right edge?  I.e., is the current state off the right of the alignment?
     * @return true if off the right edge, false if otherwise
     */
    public boolean isRightEdge() {
        return readOffset == read.getLength();
    }

    /**
     * What is our current offset in the read's bases that aligns us with the reference genome?
     *
     * @return the current read offset position.  If an edge will be == -1
     */
    public int getReadOffset() {
        return readOffset;
    }

    /**
     * What is the current offset w.r.t. the alignment state that aligns us to the readOffset?
     *
     * @return the current offset from the alignment start on the genome.  If this state is
     * at the left edge the result will be -1;
     */
    public int getGenomeOffset() {
        return genomeOffset;
    }

    /**
     * Get the position (1-based as standard) of the current alignment on the genome w.r.t. the read's alignment start
     * @return the position on the genome of the current state in absolute coordinates
     */
    public int getGenomePosition() {
        return read.getStart() + getGenomeOffset();
    }

    /**
     * Gets #getGenomePosition but as a 1 bp GenomeLoc
     * @return a non-null genome location with start position of getGenomePosition
     */
    public SimpleInterval getLocation() {
        // TODO -- may return wonky results if on an edge (could be 0 or could be beyond genome location)
        return new SimpleInterval(read.getContig(), getGenomePosition(), getGenomePosition());
    }

    /**
     * Get the cigar element we're currently aligning with.
     *
     * For example, if the cigar string is 2M2D2M and we're in the second step of the
     * first 2M, then this function returns the element 2M.  After calling stepForwardOnGenome
     * this function would return 2D.
     *
     * @return the cigar element, or null if we're the left edge
     */
    public CigarElement getCurrentCigarElement() {
        return currentElement;
    }

    /**
     * Get the offset of the current cigar element among all cigar elements in the read
     *
     * Suppose our read's cigar is 1M2D3M, and we're at the first 1M.  This would
     * return 0.  Stepping forward puts us in the 2D, so our offset is 1.  Another
     * step forward would result in a 1 again (we're in the second position of the 2D).
     * Finally, one more step forward brings us to 2 (for the 3M element)
     *
     * @return the offset of the current cigar element in the reads's cigar.  Will return -1 for
     * when the state is on the left edge, and be == the number of cigar elements in the
     * read when we're past the last position on the genome
     */
    public int getCurrentCigarElementOffset() {
        return currentCigarElementOffset;
    }

    /**
     * Get the offset of the current state into the current cigar element
     *
     * That is, suppose we have a read with cigar 2M3D4M, and we're right at
     * the second M position.  offsetIntoCurrentCigarElement would be 1, as
     * it's two elements into the 2M cigar.  Now stepping forward we'd be
     * in cigar element 3D, and our offsetIntoCurrentCigarElement would be 0.
     *
     * @return the offset (from 0) of the current state in the current cigar element.
     *  Will be 0 on the right edge, and -1 on the left.
     */
    public int getOffsetIntoCurrentCigarElement() {
        return offsetIntoCurrentCigarElement;
    }

    /**
     * Convenience accessor of the CigarOperator of the current cigar element
     *
     * Robust to the case where we're on the edge, and currentElement is null, in which
     * case this function returns null as well
     *
     * @return null if this is an edge state
     */
    public CigarOperator getCigarOperator() {
        return currentElement == null ? null : currentElement.getOperator();
    }

    @Override
    public String toString() {
        return String.format("%s ro=%d go=%d cec=%d %s", read.getName(), readOffset, genomeOffset, offsetIntoCurrentCigarElement, currentElement);
    }

    // -----------------------------------------------------------------------------------------------
    //
    // Code for setting up prev / next states
    //
    // -----------------------------------------------------------------------------------------------

    /**
     * Step the state machine forward one unit
     *
     * Takes the current state of this machine, and advances the state until the next on-genome
     * cigar element (M, X, =, D) is encountered, at which point this function returns with the
     * cigar operator of the current element.
     *
     * Assumes that the AlignmentStateMachine is in the left edge state at the start, so that
     * stepForwardOnGenome() can be called to move the machine to the first alignment position.  That
     * is, the normal use of this code is:
     *
     * AlignmentStateMachine machine = new AlignmentStateMachine(read)
     * machine.stepForwardOnGenome()
     * // now the machine is at the first position on the genome
     *
     * When stepForwardOnGenome() advances off the right edge of the read, the state machine is
     * left in a state such that isRightEdge() returns true and returns null, indicating the
     * the machine cannot advance further.  The machine may explode, though this is not contracted,
     * if stepForwardOnGenome() is called after a previous call returned null.
     *
     * @return the operator of the cigar element that machine stopped at, null if we advanced off the end of the read
     */
    public CigarOperator stepForwardOnGenome() {
        // loop until we either find a cigar element step that moves us one base on the genome, or we run
        // out of cigar elements
        while ( true ) {
            // we enter this method with readOffset = index of the last processed base on the read
            // (-1 if we did not process a single base yet); this can be last matching base,
            // or last base of an insertion
            if (currentElement == null || (offsetIntoCurrentCigarElement + 1) >= currentElement.getLength()) {
                currentCigarElementOffset++;
                if (currentCigarElementOffset < nCigarElements) {
                    currentElement = cigar.getCigarElement(currentCigarElementOffset);
                    offsetIntoCurrentCigarElement = -1;
                    // next line: guards against cigar elements of length 0; when new cigar element is retrieved,
                    // we reenter in order to re-check offsetIntoCurrentCigarElement against currentElement's length
                    continue;
                } else {
                    if (currentElement != null && currentElement.getOperator() == CigarOperator.D)
                        throw new UserException.MalformedRead(read, "read ends with deletion. Cigar: " + read.getCigar().toString() + ". Although the SAM spec technically permits such reads, this is often indicative of malformed files.");

                    // we're done, so set the offset of the cigar to 0 for cleanliness, as well as the current element
                    offsetIntoCurrentCigarElement = 0;
                    readOffset = read.getLength();
                    currentElement = null;

                    // Reads that contain indels model the genomeOffset as the following base in the reference.  Because
                    // we fall into this else block only when indels end the read, increment genomeOffset  such that the
                    // current offset of this read is the next ref base after the end of the indel.  This position will
                    // model a point on the reference somewhere after the end of the read.
                    genomeOffset++; // extended events need that. Logically, it's legal to advance the genomic offset here:

                    // we do step forward on the ref, and by returning null we also indicate that we are past the read end.
                    return null;
                }
            }

            offsetIntoCurrentCigarElement++;
            boolean done = false;
            switch (currentElement.getOperator()) {
                case H: // ignore hard clips
                case P: // ignore pads
                    offsetIntoCurrentCigarElement = currentElement.getLength();
                    break;
                case I: // insertion w.r.t. the reference
                case S: // soft clip
                    offsetIntoCurrentCigarElement = currentElement.getLength();
                    readOffset += currentElement.getLength();
                    break;
                case D: // deletion w.r.t. the reference
                    if (readOffset < 0)             // we don't want reads starting with deletion, this is a malformed cigar string
                        throw new UserException.MalformedRead(read, "read starts with deletion. Cigar: " + read.getCigar().toString() + ". Although the SAM spec technically permits such reads, this is often indicative of malformed files.");
                    // should be the same as N case
                    genomeOffset++;
                    done = true;
                    break;
                case N: // reference skip (looks and gets processed just like a "deletion", just different logical meaning)
                    genomeOffset++;
                    done = true;
                    break;
                case M:
                case EQ:
                case X:
                    readOffset++;
                    genomeOffset++;
                    done = true;
                    break;
                default:
                    throw new IllegalStateException("Case statement didn't deal with cigar op: " + currentElement.getOperator());
            }

            if ( done )
                return currentElement.getOperator();
        }
    }

    /**
     * Create a new PileupElement based on the current state of this element
     *
     * Must not be a left or right edge
     *
     * @return a pileup element
     */
    public final PileupElement makePileupElement() {
        Utils.validate(!(isLeftEdge() || isRightEdge()), "Cannot make a pileup element from an edge alignment state");
        return new PileupElement(read,
                getReadOffset(),
                getCurrentCigarElement(),
                getCurrentCigarElementOffset(),
                getOffsetIntoCurrentCigarElement());
    }
}

