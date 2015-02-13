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

package org.broadinstitute.hellbender.utils.clipping;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.ClippingRepresentation;
import org.broadinstitute.hellbender.utils.sam.CigarUtils;
import org.broadinstitute.hellbender.utils.sam.ReadUtils;

import java.util.ArrayList;
import java.util.List;

import static org.broadinstitute.hellbender.utils.sam.ReadUtils.*;

/**
 * A comprehensive clipping tool.
 *
 * General Contract:
 *  - All clipping operations return a new read with the clipped bases requested, it never modifies the original read.
 *  - If a read is fully clipped, return an empty SAMRecord, never null.
 *  - When hard clipping, add cigar operator H for every *reference base* removed (i.e. Matches, SoftClips and Deletions, but *not* insertions). See Hard Clipping notes for details.
 *
 *
 * There are several types of clipping to use:
 *
 * Write N's:
 *   Change the bases to N's in the desired region. This can be applied anywhere in the read.
 *
 * Write Q0's:
 *   Change the quality of the bases in the desired region to Q0. This can be applied anywhere in the read.
 *
 * Write both N's and Q0's:
 *   Same as the two independent operations, put together.
 *
 * Soft Clipping:
 *   Do not change the read, just mark the reads as soft clipped in the Cigar String
 *   and adjust the alignment start and end of the read.
 *
 * Hard Clipping:
 *   Creates a new read without the hard clipped bases (and base qualities). The cigar string
 *   will be updated with the cigar operator H for every reference base removed (i.e. Matches,
 *   Soft clipped bases and deletions, but *not* insertions). This contract with the cigar
 *   is necessary to allow read.getUnclippedStart() / End() to recover the original alignment
 *   of the read (before clipping).
 *
 */
public class ReadClipper {
    final SAMRecord read;
    boolean wasClipped;
    List<org.broadinstitute.hellbender.utils.ClippingOp> ops = null;

    /**
     * Initializes a ReadClipper object.
     *
     * You can set up your clipping operations using the addOp method. When you're ready to
     * generate a new read with all the clipping operations, use clipRead().
     *
     * Note: Use this if you want to set up multiple operations on the read using the ClippingOp
     * class. If you just want to apply one of the typical modes of clipping, use the static
     * clipping functions available in this class instead.
     *
     * @param read the read to clip
     */
    public ReadClipper(final SAMRecord read) {
        this.read = read;
        this.wasClipped = false;
    }

    /**
     * Add clipping operation to the read.
     *
     * You can add as many operations as necessary to this read before clipping. Beware that the
     * order in which you add these operations matter. For example, if you hard clip the beginning
     * of a read first then try to hard clip the end, the indices will have changed. Make sure you
     * know what you're doing, otherwise just use the static functions below that take care of the
     * ordering for you.
     *
     * Note: You only choose the clipping mode when you use clipRead()
     *
     * @param op a ClippingOp object describing the area you want to clip.
     */
    public void addOp(org.broadinstitute.hellbender.utils.ClippingOp op) {
        if (ops == null) ops = new ArrayList<org.broadinstitute.hellbender.utils.ClippingOp>();
        ops.add(op);
    }

    /**
     * Check the list of operations set up for this read.
     *
     * @return a list of the operations set up for this read.
     */
    public List<org.broadinstitute.hellbender.utils.ClippingOp> getOps() {
        return ops;
    }

    /**
     * Check whether or not this read has been clipped.
     * @return true if this read has produced a clipped read, false otherwise.
     */
    public boolean wasClipped() {
        return wasClipped;
    }

    /**
     * The original read.
     *
     * @return  returns the read to be clipped (original)
     */
    public SAMRecord getRead() {
        return read;
    }

    /**
     * Clips a read according to ops and the chosen algorithm.
     *
     * @param algorithm What mode of clipping do you want to apply for the stacked operations.
     * @return the read with the clipping applied.
     */
    public SAMRecord clipRead(org.broadinstitute.hellbender.utils.ClippingRepresentation algorithm) {
        if (ops == null)
            return getRead();

        SAMRecord clippedRead = read;
        for (org.broadinstitute.hellbender.utils.ClippingOp op : getOps()) {
            final int readLength = clippedRead.getReadLength();
            //check if the clipped read can still be clipped in the range requested
            if (op.start < readLength) {
                org.broadinstitute.hellbender.utils.ClippingOp fixedOperation = op;
                if (op.stop >= readLength)
                    fixedOperation = new org.broadinstitute.hellbender.utils.ClippingOp(op.start, readLength - 1);

                clippedRead = fixedOperation.apply(algorithm, clippedRead);
            }
        }
        wasClipped = true;
        ops.clear();
        if ( isEmpty(clippedRead) )
            return ReadUtils.emptyRead(clippedRead);
        return clippedRead;
    }


    /**
     * Hard clips the left tail of a read up to (and including) refStop using reference
     * coordinates.
     *
     * @param refStop the last base to be hard clipped in the left tail of the read.
     * @return a new read, without the left tail.
     */
    private SAMRecord hardClipByReferenceCoordinatesLeftTail(int refStop) {
        return hardClipByReferenceCoordinates(-1, refStop);
    }
    public static SAMRecord hardClipByReferenceCoordinatesLeftTail(SAMRecord read, int refStop) {
        return (new ReadClipper(read)).hardClipByReferenceCoordinates(-1, refStop);
    }



    /**
     * Hard clips the right tail of a read starting at (and including) refStart using reference
     * coordinates.
     *
     * @param refStart refStop the first base to be hard clipped in the right tail of the read.
     * @return a new read, without the right tail.
     */
    private SAMRecord hardClipByReferenceCoordinatesRightTail(int refStart) {
        return hardClipByReferenceCoordinates(refStart, -1);
    }
    public static SAMRecord hardClipByReferenceCoordinatesRightTail(SAMRecord read, int refStart) {
        return (new ReadClipper(read)).hardClipByReferenceCoordinates(refStart, -1);
    }

    /**
     * Hard clips a read using read coordinates.
     *
     * @param start the first base to clip (inclusive)
     * @param stop the last base to clip (inclusive)
     * @return a new read, without the clipped bases
     */
    private SAMRecord hardClipByReadCoordinates(int start, int stop) {
        if (isEmpty(read) || (start == 0 && stop == read.getReadLength() - 1))
            return ReadUtils.emptyRead(read);

        this.addOp(new org.broadinstitute.hellbender.utils.ClippingOp(start, stop));
        return clipRead(org.broadinstitute.hellbender.utils.ClippingRepresentation.HARDCLIP_BASES);
    }
    public static SAMRecord hardClipByReadCoordinates(SAMRecord read, int start, int stop) {
        return (new ReadClipper(read)).hardClipByReadCoordinates(start, stop);
    }


    /**
     * Hard clips both tails of a read.
     *   Left tail goes from the beginning to the 'left' coordinate (inclusive)
     *   Right tail goes from the 'right' coordinate (inclusive) until the end of the read
     *
     * @param left the coordinate of the last base to be clipped in the left tail (inclusive)
     * @param right the coordinate of the first base to be clipped in the right tail (inclusive)
     * @return a new read, without the clipped bases
     */
    private SAMRecord hardClipBothEndsByReferenceCoordinates(int left, int right) {
        if (isEmpty(read) || left == right)
            return ReadUtils.emptyRead(read);
        SAMRecord leftTailRead = hardClipByReferenceCoordinates(right, -1);

        // after clipping one tail, it is possible that the consequent hard clipping of adjacent deletions
        // make the left cut index no longer part of the read. In that case, clip the read entirely.
        if (left > leftTailRead.getAlignmentEnd())
            return ReadUtils.emptyRead(read);

        ReadClipper clipper = new ReadClipper(leftTailRead);
        return clipper.hardClipByReferenceCoordinatesLeftTail(left);
    }
    public static SAMRecord hardClipBothEndsByReferenceCoordinates(SAMRecord read, int left, int right) {
        return (new ReadClipper(read)).hardClipBothEndsByReferenceCoordinates(left, right);
    }


    /**
     * Clips any contiguous tail (left, right or both) with base quality lower than lowQual using the desired algorithm.
     *
     * This function will look for low quality tails and hard clip them away. A low quality tail
     * ends when a base has base quality greater than lowQual.
     *
     * @param algorithm the algorithm to use (HardClip, SoftClip, Write N's,...)
     * @param lowQual every base quality lower than or equal to this in the tail of the read will be hard clipped
     * @return a new read without low quality tails
     */
    private SAMRecord clipLowQualEnds(org.broadinstitute.hellbender.utils.ClippingRepresentation algorithm, byte lowQual) {
        if (isEmpty(read))
            return read;

        final byte [] quals = read.getBaseQualities();
        final int readLength = read.getReadLength();
        int leftClipIndex = 0;
        int rightClipIndex = readLength - 1;

        // check how far we can clip both sides
        while (rightClipIndex >= 0 && quals[rightClipIndex] <= lowQual) rightClipIndex--;
        while (leftClipIndex < readLength && quals[leftClipIndex] <= lowQual) leftClipIndex++;

        // if the entire read should be clipped, then return an empty read.
        if (leftClipIndex > rightClipIndex)
            return ReadUtils.emptyRead(read);

        if (rightClipIndex < readLength - 1) {
            this.addOp(new org.broadinstitute.hellbender.utils.ClippingOp(rightClipIndex + 1, readLength - 1));
        }
        if (leftClipIndex > 0 ) {
            this.addOp(new org.broadinstitute.hellbender.utils.ClippingOp(0, leftClipIndex - 1));
        }
        return this.clipRead(algorithm);
    }

    private SAMRecord hardClipLowQualEnds(byte lowQual) {
        return this.clipLowQualEnds(org.broadinstitute.hellbender.utils.ClippingRepresentation.HARDCLIP_BASES, lowQual);
    }
    public static SAMRecord hardClipLowQualEnds(SAMRecord read, byte lowQual) {
        return (new ReadClipper(read)).hardClipLowQualEnds(lowQual);
    }
    public static SAMRecord clipLowQualEnds(SAMRecord read, byte lowQual, org.broadinstitute.hellbender.utils.ClippingRepresentation algorithm) {
        return (new ReadClipper(read)).clipLowQualEnds(algorithm, lowQual);
    }


    /**
     * Will hard clip every soft clipped bases in the read.
     *
     * @return a new read without the soft clipped bases
     */
    private SAMRecord hardClipSoftClippedBases () {
        if (isEmpty(read))
            return read;

        int readIndex = 0;
        int cutLeft = -1;            // first position to hard clip (inclusive)
        int cutRight = -1;           // first position to hard clip (inclusive)
        boolean rightTail = false;   // trigger to stop clipping the left tail and start cutting the right tail

        for (CigarElement cigarElement : read.getCigar().getCigarElements()) {
            if (cigarElement.getOperator() == CigarOperator.SOFT_CLIP) {
                if (rightTail) {
                    cutRight = readIndex;
                }
                else {
                    cutLeft = readIndex + cigarElement.getLength() - 1;
                }
            }
            else if (cigarElement.getOperator() != CigarOperator.HARD_CLIP)
                rightTail = true;

            if (cigarElement.getOperator().consumesReadBases())
                readIndex += cigarElement.getLength();
        }

        // It is extremely important that we cut the end first otherwise the read coordinates change.
        if (cutRight >= 0)
            this.addOp(new org.broadinstitute.hellbender.utils.ClippingOp(cutRight, read.getReadLength() - 1));
        if (cutLeft >= 0)
            this.addOp(new org.broadinstitute.hellbender.utils.ClippingOp(0, cutLeft));

        return clipRead(org.broadinstitute.hellbender.utils.ClippingRepresentation.HARDCLIP_BASES);
    }
    public static SAMRecord hardClipSoftClippedBases (SAMRecord read) {
        return (new ReadClipper(read)).hardClipSoftClippedBases();
    }


    /**
     * Hard clip the read to the variable region (from refStart to refStop)
     *
     * @param read     the read to be clipped
     * @param refStart the beginning of the variant region (inclusive)
     * @param refStop  the end of the variant region (inclusive)
     * @return the read hard clipped to the variant region
     */
    public static SAMRecord hardClipToRegion( final SAMRecord read, final int refStart, final int refStop ) {
        final int start = read.getAlignmentStart();
        final int stop = read.getAlignmentEnd();
        return hardClipToRegion(read, refStart, refStop,start,stop);
    }

    /**
     * Hard clip the read to the variable region (from refStart to refStop) processing also the clipped bases
     *
     * @param read     the read to be clipped
     * @param refStart the beginning of the variant region (inclusive)
     * @param refStop  the end of the variant region (inclusive)
     * @return the read hard clipped to the variant region
     */
    public static SAMRecord hardClipToRegionIncludingClippedBases( final SAMRecord read, final int refStart, final int refStop ) {
        final int start = read.getUnclippedStart();
        final int stop = start + CigarUtils.countRefBasesBasedOnCigar(read, 0, read.getCigarLength()) - 1;
        return hardClipToRegion(read, refStart, refStop,start,stop);
    }

    private static SAMRecord hardClipToRegion( final SAMRecord read, final int refStart, final int refStop, final int alignmentStart, final int alignmentStop){
        // check if the read is contained in region
        if (alignmentStart <= refStop && alignmentStop >= refStart) {
            if (alignmentStart < refStart && alignmentStop > refStop)
                return hardClipBothEndsByReferenceCoordinates(read, refStart - 1, refStop + 1);
            else if (alignmentStart < refStart)
                return hardClipByReferenceCoordinatesLeftTail(read, refStart - 1);
            else if (alignmentStop > refStop)
                return hardClipByReferenceCoordinatesRightTail(read, refStop + 1);
            return read;
        } else
            return ReadUtils.emptyRead(read);

    }

    public static List<SAMRecord> hardClipToRegion( final List<SAMRecord> reads, final int refStart, final int refStop ) {
        final List<SAMRecord> returnList = new ArrayList<SAMRecord>( reads.size() );
        for( final SAMRecord read : reads ) {
            final SAMRecord clippedRead = hardClipToRegion( read, refStart, refStop );
            if( !isEmpty(clippedRead) ) {
                returnList.add( clippedRead );
            }
        }
        return returnList;
    }

    /**
     * Checks if a read contains adaptor sequences. If it does, hard clips them out.
     *
     * Note: To see how a read is checked for adaptor sequence see ReadUtils.getAdaptorBoundary()
     *
     * @return a new read without adaptor sequence
     */
    private SAMRecord hardClipAdaptorSequence () {
        final int adaptorBoundary = getAdaptorBoundary(read);

        if (adaptorBoundary == CANNOT_COMPUTE_ADAPTOR_BOUNDARY || !isInsideRead(read, adaptorBoundary))
            return read;

        return read.getReadNegativeStrandFlag() ? hardClipByReferenceCoordinatesLeftTail(adaptorBoundary) : hardClipByReferenceCoordinatesRightTail(adaptorBoundary);
    }
    public static SAMRecord hardClipAdaptorSequence (SAMRecord read) {
        return (new ReadClipper(read)).hardClipAdaptorSequence();
    }


    /**
     * Hard clips any leading insertions in the read. Only looks at the beginning of the read, not the end.
     *
     * @return a new read without leading insertions
     */
    private SAMRecord hardClipLeadingInsertions() {
        if (isEmpty(read))
            return read;

        for(CigarElement cigarElement : read.getCigar().getCigarElements()) {
            if (cigarElement.getOperator() != CigarOperator.HARD_CLIP && cigarElement.getOperator() != CigarOperator.SOFT_CLIP &&
                    cigarElement.getOperator() != CigarOperator.INSERTION)
                break;

            else if (cigarElement.getOperator() == CigarOperator.INSERTION)
                this.addOp(new org.broadinstitute.hellbender.utils.ClippingOp(0, cigarElement.getLength() - 1));

        }
        return clipRead(org.broadinstitute.hellbender.utils.ClippingRepresentation.HARDCLIP_BASES);
    }
    public static SAMRecord hardClipLeadingInsertions(SAMRecord read) {
        return (new ReadClipper(read)).hardClipLeadingInsertions();
    }


    /**
     * Turns soft clipped bases into matches
     * @return a new read with every soft clip turned into a match
     */
    private SAMRecord revertSoftClippedBases() {
        if (isEmpty(read))
            return read;

        this.addOp(new org.broadinstitute.hellbender.utils.ClippingOp(0, 0));
        return this.clipRead(org.broadinstitute.hellbender.utils.ClippingRepresentation.REVERT_SOFTCLIPPED_BASES);
    }

    /**
     * Reverts ALL soft-clipped bases
     *
     * @param read the read
     * @return the read with all soft-clipped bases turned into matches
     */
    public static SAMRecord revertSoftClippedBases(SAMRecord read) {
        return (new ReadClipper(read)).revertSoftClippedBases();
    }

    /**
     * Reverts only soft clipped bases with quality score greater than or equal to minQual
     *
     * todo -- Note: Will write a temporary field with the number of soft clips that were undone on each side (left: 'SL', right: 'SR') -- THIS HAS BEEN REMOVED TEMPORARILY SHOULD HAPPEN INSIDE THE CLIPPING ROUTINE!
     *
     * @param read    the read
     * @param minQual the mininum base quality score to revert the base (inclusive)
     * @return a new read with high quality soft clips reverted
     */
    public static SAMRecord revertSoftClippedBases(SAMRecord read, byte minQual) {
        return revertSoftClippedBases(hardClipLowQualitySoftClips(read, minQual));
    }

    /**
     * Hard clips away soft clipped bases that are below the given quality threshold
     *
     * @param read    the read
     * @param minQual the minimum base quality score to revert the base (inclusive)
     * @return a new read without low quality soft clipped bases
     */
    public static SAMRecord hardClipLowQualitySoftClips(SAMRecord read, byte minQual) {
        int nLeadingSoftClips = read.getAlignmentStart() - getSoftStart(read);
        if (isEmpty(read) || nLeadingSoftClips > read.getReadLength())
            return ReadUtils.emptyRead(read);

        byte [] quals = read.getBaseQualities();
        int left = -1;

        if (nLeadingSoftClips > 0) {
            for (int i = nLeadingSoftClips - 1; i >= 0; i--) {
                if (quals[i] >= minQual)
                    left = i;
                else
                    break;
            }
        }

        int right = -1;
        int nTailingSoftClips = getSoftEnd(read) - read.getAlignmentEnd();
        if (nTailingSoftClips > 0) {
            for (int i = read.getReadLength() - nTailingSoftClips; i < read.getReadLength() ; i++) {
                if (quals[i] >= minQual)
                    right = i;
                else
                    break;
            }
        }

        SAMRecord clippedRead = read;
        if (right >= 0 && right + 1 < clippedRead.getReadLength())                                                      // only clip if there are softclipped bases (right >= 0) and the first high quality soft clip is not the last base (right+1 < readlength)
            clippedRead = hardClipByReadCoordinates(clippedRead, right+1, clippedRead.getReadLength()-1);           // first we hard clip the low quality soft clips on the right tail
        if (left >= 0 && left - 1 > 0)                                                                                  // only clip if there are softclipped bases (left >= 0) and the first high quality soft clip is not the last base (left-1 > 0)
            clippedRead = hardClipByReadCoordinates(clippedRead, 0, left-1);                                        // then we hard clip the low quality soft clips on the left tail

        return clippedRead;
    }

    /**
     * Generic functionality to hard clip a read, used internally by hardClipByReferenceCoordinatesLeftTail
     * and hardClipByReferenceCoordinatesRightTail. Should not be used directly.
     *
     * Note, it REQUIRES you to give the directionality of your hard clip (i.e. whether you're clipping the
     * left of right tail) by specifying either refStart < 0 or refStop < 0.
     *
     * @param refStart  first base to clip (inclusive)
     * @param refStop last base to clip (inclusive)
     * @return a new read, without the clipped bases
     */
    protected SAMRecord hardClipByReferenceCoordinates(int refStart, int refStop) {
        if (isEmpty(read))
            return read;

        int start;
        int stop;

        // Determine the read coordinate to start and stop hard clipping
        if (refStart < 0) {
            if (refStop < 0)
                throw new GATKException("Only one of refStart or refStop must be < 0, not both (" + refStart + ", " + refStop + ")");
            start = 0;
            stop = ReadUtils.getReadCoordinateForReferenceCoordinate(read, refStop, ReadUtils.ClippingTail.LEFT_TAIL);
        }
        else {
            if (refStop >= 0)
                throw new GATKException("Either refStart or refStop must be < 0 (" + refStart + ", " + refStop + ")");
            start = ReadUtils.getReadCoordinateForReferenceCoordinate(read, refStart, ReadUtils.ClippingTail.RIGHT_TAIL);
            stop = read.getReadLength() - 1;
        }

        if (start < 0 || stop > read.getReadLength() - 1)
            throw new GATKException("Trying to clip before the start or after the end of a read");

        if ( start > stop )
            throw new GATKException(String.format("START (%d) > (%d) STOP -- this should never happen, please check read: %s (CIGAR: %s)", start, stop, read, read.getCigarString()));

        if ( start > 0 && stop < read.getReadLength() - 1)
            throw new GATKException(String.format("Trying to clip the middle of the read: start %d, stop %d, cigar: %s", start, stop, read.getCigarString()));

        this.addOp(new org.broadinstitute.hellbender.utils.ClippingOp(start, stop));
        SAMRecord clippedRead = clipRead(ClippingRepresentation.HARDCLIP_BASES);
        this.ops = null;
        return clippedRead;
    }
}