package org.broadinstitute.hellbender.utils.clipping;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;

/**
 * Represents a clip on a read.  It has a type (see the enum) along with a start and stop in the bases
 * of the read, plus an option extraInfo (useful for carrying info where needed).
 * <p/>
 * Also holds the critical apply function that actually execute the clipping operation on a provided read,
 * according to the wishes of the supplied ClippingAlgorithm enum.
 */
public final class ClippingOp {
    public final int start, stop; // inclusive

    public ClippingOp(final int start, final int stop) {
        this.start = start;
        this.stop = stop;
    }


    public int getLength() {
        return stop - start + 1;
    }

    /**
     * Clips the bases in read according to this operation's start and stop.  Uses the clipping
     * representation used is the one provided by algorithm argument.
     *
     * @param algorithm    clipping algorithm to use
     * @param originalRead the read to be clipped
     */
    public GATKRead apply(final ClippingRepresentation algorithm, final GATKRead originalRead) {
        return apply(algorithm, originalRead, true);
    }

    // Warning, this method is intended for testing purposes
    GATKRead apply(final ClippingRepresentation algorithm, final GATKRead originalRead, final boolean runAsserts) {
        switch (algorithm) {
            // important note:
            //   it's not safe to call read.getBases()[i] = 'N' or read.getBaseQualities()[i] = 0
            //   because you're not guaranteed to get a pointer to the actual array of bytes in the GATKRead
            case WRITE_NS: {
                final GATKRead readCopied = originalRead.copy();
                applyWriteNs(readCopied);
                return readCopied;
            }

            case WRITE_Q0S: {
                final GATKRead readCopied = originalRead.copy();
                applyWriteQ0s(readCopied);
                return readCopied;
            }

            case WRITE_NS_Q0S: {
                final GATKRead readCopied = originalRead.copy();
                applyWriteNs(readCopied);
                applyWriteQ0s(readCopied);
                return readCopied;
            }

            case HARDCLIP_BASES: {
                //Note: passing the original read here because the read is copied right away inside the method
                return applyHardClipBases(originalRead, start, stop);
            }

            case SOFTCLIP_BASES: {
                return applySoftClipBases(originalRead.copy(), runAsserts);
            }

            case REVERT_SOFTCLIPPED_BASES: {
                return applyRevertSoftClippedBases(originalRead.copy());
            }

            default: {
                throw new IllegalStateException("Unexpected Clipping operator type " + algorithm);
            }
        }
    }

    private GATKRead applySoftClipBases(final GATKRead readCopied, final boolean runAsserts) {
        if (readCopied.isUnmapped()) {
            // we can't process unmapped reads
            throw new UserException("Read Clipper cannot soft clip unmapped reads");
        }

        int myStop = stop;
        if ((stop + 1 - start) == readCopied.getLength()) {
            // BAM representation issue -- we can't SOFTCLIP away all bases in a read, just leave it alone
            myStop--; // just decrement stop
        }

        if (start > 0 && myStop != readCopied.getLength() - 1) {
            throw new GATKException(String.format("Cannot apply soft clipping operator to the middle of a read: %s to be clipped at %d-%d", readCopied.getName(), start, myStop));
        }

        final Cigar oldCigar = readCopied.getCigar();

        int scLeft = 0;
        int scRight = readCopied.getLength();
        if (start == 0) {
            scLeft = myStop + 1;
        } else {
            scRight = start;
        }

        final Cigar newCigar = softClip(oldCigar, scLeft, scRight, runAsserts);
        readCopied.setCigar(newCigar);

        final int newClippedStart = getNewAlignmentStartOffset(newCigar, oldCigar);
        final int newStart = readCopied.getStart() + newClippedStart;
        readCopied.setPosition(readCopied.getContig(), newStart);
        return readCopied;
    }

    private void applyWriteQ0s(final GATKRead readCopied) {
        final byte[] newQuals = readCopied.getBaseQualities(); //this makes a copy so we can modify in place
        overwriteFromStartToStop(newQuals, (byte)0);
        readCopied.setBaseQualities(newQuals);
    }

    private void applyWriteNs(final GATKRead readCopied) {
        final byte[] newBases = readCopied.getBases();       //this makes a copy so we can modify in place
        overwriteFromStartToStop(newBases, (byte)'N');
        readCopied.setBases(newBases);
    }

    private void overwriteFromStartToStop(final byte[] arr, final byte newVal){
        Arrays.fill(arr, start, Math.min(arr.length, stop+1), newVal);
    }

    private GATKRead applyRevertSoftClippedBases(final GATKRead read) {
        GATKRead unclipped = read.copy();

        final Cigar unclippedCigar = new Cigar();
        int matchesCount = 0;
        for (final CigarElement element : read.getCigarElements()) {
            if (element.getOperator() == CigarOperator.SOFT_CLIP || element.getOperator() == CigarOperator.MATCH_OR_MISMATCH) {
                matchesCount += element.getLength();
            } else if (matchesCount > 0) {
                unclippedCigar.add(new CigarElement(matchesCount, CigarOperator.MATCH_OR_MISMATCH));
                matchesCount = 0;
                unclippedCigar.add(element);
            } else {
                unclippedCigar.add(element);
            }
        }
        if (matchesCount > 0) {
            unclippedCigar.add(new CigarElement(matchesCount, CigarOperator.MATCH_OR_MISMATCH));
        }

        unclipped.setCigar(unclippedCigar);
        final int newStart = read.getStart() + calculateAlignmentStartShift(read.getCigar(), unclippedCigar);

        if ( newStart <= 0 ) {
            // if the start of the unclipped read occurs before the contig,
            // we must hard clip away the bases since we cannot represent reads with
            // negative or 0 alignment start values in the SAMRecord (e.g., 0 means unaligned)

            // We cannot set the read to temporarily have a negative start position, as our Read
            // interface will not allow it. Instead, since we know that the final start position will
            // be 1 after the hard clip operation, set it to 1 explicitly. We have to set it twice:
            // once before the hard clip (to reset the alignment stop / read length in read implementations
            // that cache these values, such as SAMRecord), and again after the hard clip.
            unclipped.setPosition(unclipped.getContig(), 1);
            unclipped = applyHardClipBases(unclipped, 0, - newStart);
            
            // Reset the position to 1 again only if we didn't end up with an empty, unmapped read after hard clipping.
            // See https://github.com/broadinstitute/gatk/issues/3845
            if ( ! unclipped.isUnmapped() ) {
                unclipped.setPosition(unclipped.getContig(), 1);
            }
            
            return unclipped;
        } else {
            unclipped.setPosition(unclipped.getContig(), newStart);
            return unclipped;
        }
    }

    /**
     * Given two cigar strings corresponding to read before and after soft-clipping, returns an integer
     * corresponding to the number of reference bases that the new string must be offset by in order to have the
     * correct start according to the reference.
     *
     * @param clippedCigar the new cigar string after clipping
     * @param oldCigar the cigar string of the read before it was soft clipped
     */
    @VisibleForTesting
    static int getNewAlignmentStartOffset(final Cigar clippedCigar, final Cigar oldCigar) {
        int readBasesBeforeReference = 0; // The number of read bases consumed on the new cigar before reference bases are consumed

        int basesBeforeReferenceOld = 0; // The number of read bases consumed on the old cigar before reference bases are consumed
        int curReadCounter = 0; // A measure of the reference offset between the oldCigar and the clippedCigar

        for (final CigarElement e : clippedCigar.getCigarElements()) {
            if (!e.getOperator().consumesReferenceBases()) {
                if (e.getOperator().consumesReadBases()) {
                    readBasesBeforeReference += e.getLength();
                }
            } else {
                if (!e.getOperator().consumesReadBases()) {
                    basesBeforeReferenceOld -= e.getLength(); // Accounting for any D or N cigar operators at the front of the string
                } else {
                    break;
                }
            }
        }


        for (final CigarElement e : oldCigar.getCigarElements()) {
            int curRefLength = e.getLength();
            int curReadLength = e.getLength();
            if (!e.getOperator().consumesReadBases()) {
                curReadLength = 0;
            }

            boolean truncated = false;
            if (curReadCounter + curReadLength > readBasesBeforeReference) {
                curReadLength = readBasesBeforeReference - curReadCounter;
                curRefLength = readBasesBeforeReference - curReadCounter;
                truncated = true;
            }

            if (!e.getOperator().consumesReferenceBases()) {
                curRefLength = 0;
            }

            curReadCounter += curReadLength;
            basesBeforeReferenceOld += curRefLength;

            if (curReadCounter > readBasesBeforeReference || truncated) {
                break;
            }
        }

        return Math.abs(basesBeforeReferenceOld); // if oldNum is negative it means some of the preceding N/Ds were trimmed but not all so we take absolute value
    }

    /**
     * Given a cigar string, soft clip up to startClipEnd and soft clip starting at endClipBegin
     */
    private Cigar softClip(final Cigar __cigar, final int __startClipEnd, final int __endClipBegin, final boolean runAsserts) {
        if (__endClipBegin <= __startClipEnd) {
            //whole thing should be soft clipped
            int cigarLength = 0;
            for (final CigarElement e : __cigar.getCigarElements()) {
                cigarLength += e.getLength();
            }

            final Cigar newCigar = new Cigar();
            newCigar.add(new CigarElement(cigarLength, CigarOperator.SOFT_CLIP));
            if (runAsserts) {
                assert newCigar.isValid(null, -1) == null;
            }
            return newCigar;
        }

        int curLength = 0;
        final Vector<CigarElement> newElements = new Vector<>();
        for (final CigarElement curElem : __cigar.getCigarElements()) {
            if (!curElem.getOperator().consumesReadBases()) {
                if (curElem.getOperator() == CigarOperator.HARD_CLIP || curLength > __startClipEnd && curLength < __endClipBegin) {
                    newElements.add(new CigarElement(curElem.getLength(), curElem.getOperator()));
                }
                continue;
            }

            final int s = curLength;
            final int e = curLength + curElem.getLength();
            if (e <= __startClipEnd || s >= __endClipBegin) {
                //must turn this entire thing into a clip
                newElements.add(new CigarElement(curElem.getLength(), CigarOperator.SOFT_CLIP));
            } else if (s >= __startClipEnd && e <= __endClipBegin) {
                //same thing
                newElements.add(new CigarElement(curElem.getLength(), curElem.getOperator()));
            } else {
                //we are clipping in the middle of this guy
                CigarElement newStart = null;
                CigarElement newMid = null;
                CigarElement newEnd = null;

                int midLength = curElem.getLength();
                if (s < __startClipEnd) {
                    newStart = new CigarElement(__startClipEnd - s, CigarOperator.SOFT_CLIP);
                    midLength -= newStart.getLength();
                }

                if (e > __endClipBegin) {
                    newEnd = new CigarElement(e - __endClipBegin, CigarOperator.SOFT_CLIP);
                    midLength -= newEnd.getLength();
                }
                if (runAsserts) {
                    assert midLength >= 0;
                }
                if (midLength > 0) {
                    newMid = new CigarElement(midLength, curElem.getOperator());
                }
                if (newStart != null) {
                    newElements.add(newStart);
                }
                if (newMid != null) {
                    newElements.add(newMid);
                }
                if (newEnd != null) {
                    newElements.add(newEnd);
                }
            }
            curLength += curElem.getLength();
        }

        final Vector<CigarElement> finalNewElements = new Vector<>();
        CigarElement lastElement = null;
        for (final CigarElement elem : newElements) {
            if (lastElement == null || lastElement.getOperator() != elem.getOperator()) {
                if (lastElement != null) {
                    finalNewElements.add(lastElement);
                }
                lastElement = elem;
            } else {
                lastElement = new CigarElement(lastElement.getLength() + elem.getLength(), lastElement.getOperator());
            }
        }
        if (lastElement != null) {
            finalNewElements.add(lastElement);
        }

        final Cigar newCigar = new Cigar(finalNewElements);
        if (runAsserts) {
            assert newCigar.isValid(null, -1) == null;
        }
        return newCigar;
    }

    /**
     * Hard clip bases from read, from start to stop in base coordinates
     *
     * If start == 0, then we will clip from the front of the read, otherwise we clip
     * from the right.  If start == 0 and stop == 10, this would clip out the first
     * 10 bases of the read.
     *
     * Note that this function works with reads with negative alignment starts, in order to
     * allow us to hardClip reads that have had their soft clips reverted and so might have
     * negative alignment starts
     *
     * Works properly with reduced reads and insertion/deletion base qualities
     *
     * Note: this method does not assume that the read is directly modifiable
     * and makes a copy of it.
     *
     * @param read a non-null read
     * @param start a start >= 0 and < read.length
     * @param stop a stop >= 0 and < read.length.
     * @return a cloned version of read that has been properly trimmed down (Could be an empty, unmapped read)
     */
    private GATKRead applyHardClipBases(final GATKRead read, final int start, final int stop) {
        // If the read is unmapped there is no Cigar string and neither should we create a new cigar string

        final Cigar cigar = read.getCigar();//Get the cigar once to avoid multiple calls because each makes a copy of the cigar
        final CigarShift cigarShift = read.isUnmapped() ? new CigarShift(new Cigar(), 0, 0) : hardClipCigar(cigar, start, stop);

        // the cigar may force a shift left or right (or both) in case we are left with insertions
        // starting or ending the read after applying the hard clip on start/stop.
        final int newLength = read.getLength() - (stop - start + 1) - cigarShift.shiftFromStart - cigarShift.shiftFromEnd;

        // If the new read is going to be empty, return an empty read now. This avoids initializing the new
        // read with invalid values below in certain cases (such as a negative alignment start).
        // See https://github.com/broadinstitute/gatk/issues/3466
        if ( newLength == 0 ) {
            return ReadUtils.emptyRead(read);
        }

        final byte[] newBases = new byte[newLength];
        final byte[] newQuals = new byte[newLength];
        final int copyStart = (start == 0) ? stop + 1 + cigarShift.shiftFromStart : cigarShift.shiftFromStart;

        System.arraycopy(read.getBases(), copyStart, newBases, 0, newLength);
        System.arraycopy(read.getBaseQualities(), copyStart, newQuals, 0, newLength);

        final GATKRead hardClippedRead = read.copy();

        hardClippedRead.setBaseQualities(newQuals);
        hardClippedRead.setBases(newBases);
        hardClippedRead.setCigar(cigarShift.cigar);
        if (start == 0 && !read.isUnmapped()) {
            hardClippedRead.setPosition(read.getContig(), read.getStart() + calculateAlignmentStartShift(cigar, cigarShift.cigar));
        }

        if (ReadUtils.hasBaseIndelQualities(read)) {
            final byte[] newBaseInsertionQuals = new byte[newLength];
            final byte[] newBaseDeletionQuals = new byte[newLength];
            System.arraycopy(ReadUtils.getBaseInsertionQualities(read), copyStart, newBaseInsertionQuals, 0, newLength);
            System.arraycopy(ReadUtils.getBaseDeletionQualities(read), copyStart, newBaseDeletionQuals, 0, newLength);
            ReadUtils.setInsertionBaseQualities(hardClippedRead, newBaseInsertionQuals);
            ReadUtils.setDeletionBaseQualities(hardClippedRead, newBaseDeletionQuals);
        }

        return hardClippedRead;

    }

    private CigarShift  hardClipCigar(final Cigar cigar, final int start, final int stop) {
        final Cigar newCigar = new Cigar();
        int index = 0;
        int totalHardClipCount = stop - start + 1;
        int alignmentShift = 0; // caused by hard clipping deletions

        // hard clip the beginning of the cigar string
        if (start == 0) {
            final Iterator<CigarElement> cigarElementIterator = cigar.getCigarElements().iterator();
            CigarElement cigarElement = cigarElementIterator.next();
            // Skip all leading hard clips
            while (cigarElement.getOperator() == CigarOperator.HARD_CLIP) {
                totalHardClipCount += cigarElement.getLength();
                if (cigarElementIterator.hasNext()) {
                    cigarElement = cigarElementIterator.next();
                } else {
                    throw new GATKException("Read is entirely hard-clipped, shouldn't be trying to clip it's cigar string");
                }
            }
            // keep clipping until we hit stop
            while (index <= stop) {
                int shift = 0;
                if (cigarElement.getOperator().consumesReadBases()) {
                    shift = cigarElement.getLength();
                }

                // we're still clipping or just finished perfectly
                if (index + shift == stop + 1) {
                    alignmentShift += calculateHardClippingAlignmentShift(cigarElement, cigarElement.getLength());
                    newCigar.add(new CigarElement(totalHardClipCount + alignmentShift, CigarOperator.HARD_CLIP));
                }
                // element goes beyond what we need to clip
                else if (index + shift > stop + 1) {
                    final int elementLengthAfterChopping = cigarElement.getLength() - (stop - index + 1);
                    alignmentShift += calculateHardClippingAlignmentShift(cigarElement, stop - index + 1);
                    newCigar.add(new CigarElement(totalHardClipCount + alignmentShift, CigarOperator.HARD_CLIP));
                    newCigar.add(new CigarElement(elementLengthAfterChopping, cigarElement.getOperator()));
                }
                index += shift;
                alignmentShift += calculateHardClippingAlignmentShift(cigarElement, shift);

                if (index <= stop && cigarElementIterator.hasNext()) {
                    cigarElement = cigarElementIterator.next();
                } else {
                    break;
                }
            }

            // add the remaining cigar elements
            while (cigarElementIterator.hasNext()) {
                cigarElement = cigarElementIterator.next();
                newCigar.add(new CigarElement(cigarElement.getLength(), cigarElement.getOperator()));
            }
        }

        // hard clip the end of the cigar string
        else {
            final Iterator<CigarElement> cigarElementIterator = cigar.getCigarElements().iterator();
            CigarElement cigarElement = cigarElementIterator.next();

            // Keep marching on until we find the start
            while (index < start) {
                int shift = 0;
                if (cigarElement.getOperator().consumesReadBases()) {
                    shift = cigarElement.getLength();
                }

                // we haven't gotten to the start yet, keep everything as is.
                if (index + shift < start) {
                    newCigar.add(new CigarElement(cigarElement.getLength(), cigarElement.getOperator()));
                }// element goes beyond our clip starting position
                else {
                    final int elementLengthAfterChopping = start - index;
                    alignmentShift += calculateHardClippingAlignmentShift(cigarElement, cigarElement.getLength() - (start - index));

                    // if this last element is a HARD CLIP operator, just merge it with our hard clip operator to be added later
                    if (cigarElement.getOperator() == CigarOperator.HARD_CLIP) {
                        totalHardClipCount += elementLengthAfterChopping;
                    }// otherwise, maintain what's left of this last operator
                    else {
                        newCigar.add(new CigarElement(elementLengthAfterChopping, cigarElement.getOperator()));
                    }
                }
                index += shift;
                if (index < start && cigarElementIterator.hasNext()) {
                    cigarElement = cigarElementIterator.next();
                } else {
                    break;
                }
            }

            // check if we are hard clipping indels
            while (cigarElementIterator.hasNext()) {
                cigarElement = cigarElementIterator.next();
                alignmentShift += calculateHardClippingAlignmentShift(cigarElement, cigarElement.getLength());

                // if the read had a HardClip operator in the end, combine it with the Hard Clip we are adding
                if (cigarElement.getOperator() == CigarOperator.HARD_CLIP) {
                    totalHardClipCount += cigarElement.getLength();
                }
            }
            newCigar.add(new CigarElement(totalHardClipCount + alignmentShift, CigarOperator.HARD_CLIP));
        }
        return cleanHardClippedCigar(newCigar);
    }

    /**
     * Checks if a hard clipped cigar left a read starting or ending with deletions or gap (N)
     * and cleans it up accordingly.
     *
     * @param cigar the original cigar
     * @return an object with the shifts (see CigarShift class)
     */
    private CigarShift cleanHardClippedCigar(final Cigar cigar) {
        final Cigar cleanCigar = new Cigar();
        int shiftFromStart = 0;
        int shiftFromEnd = 0;
        Stack<CigarElement> cigarStack = new Stack<>();
        final Stack<CigarElement> inverseCigarStack = new Stack<>();

        for (final CigarElement cigarElement : cigar.getCigarElements()) {
            cigarStack.push(cigarElement);
        }

        for (int i = 1; i <= 2; i++) {
            final int shift = 0;
            int totalHardClip = 0;
            boolean readHasStarted = false;
            boolean addedHardClips = false;

            while (!cigarStack.empty()) {
                final CigarElement cigarElement = cigarStack.pop();

                if (!readHasStarted &&
                        cigarElement.getOperator() != CigarOperator.DELETION &&
                        cigarElement.getOperator() != CigarOperator.SKIPPED_REGION &&
                        cigarElement.getOperator() != CigarOperator.HARD_CLIP) {
                    readHasStarted = true;
                } else if (!readHasStarted && cigarElement.getOperator() == CigarOperator.HARD_CLIP) {
                    totalHardClip += cigarElement.getLength();
                } else if (!readHasStarted && cigarElement.getOperator() == CigarOperator.DELETION) {
                    totalHardClip += cigarElement.getLength();
                } else if (!readHasStarted && cigarElement.getOperator() == CigarOperator.SKIPPED_REGION) {
                    totalHardClip += cigarElement.getLength();
                }

                if (readHasStarted) {
                    if (i == 1) {
                        if (!addedHardClips) {
                            if (totalHardClip > 0) {
                                inverseCigarStack.push(new CigarElement(totalHardClip, CigarOperator.HARD_CLIP));
                            }
                            addedHardClips = true;
                        }
                        inverseCigarStack.push(cigarElement);
                    } else {
                        if (!addedHardClips) {
                            if (totalHardClip > 0) {
                                cleanCigar.add(new CigarElement(totalHardClip, CigarOperator.HARD_CLIP));
                            }
                            addedHardClips = true;
                        }
                        cleanCigar.add(cigarElement);
                    }
                }
            }
            // first pass  (i=1) is from end to start of the cigar elements
            if (i == 1) {
                shiftFromEnd = shift;
                cigarStack = inverseCigarStack;
            }
            // second pass (i=2) is from start to end with the end already cleaned
            else {
                shiftFromStart = shift;
            }
        }
        return new CigarShift(cleanCigar, shiftFromStart, shiftFromEnd);
    }

    /**
     * Compute the offset of the first "real" position in the cigar on the genome
     *
     * This is defined as a first position after a run of Hs followed by a run of Ss
     *
     * @param cigar A non-null cigar
     * @return the offset (from 0) of the first on-genome base
     */
    private int calcHardSoftOffset(final Cigar cigar) {
        final List<CigarElement> elements = cigar.getCigarElements();

        int size = 0;
        int i = 0;
        while ( i < elements.size() && elements.get(i).getOperator() == CigarOperator.HARD_CLIP ) {
            size += elements.get(i).getLength();
            i++;
        }
        while ( i < elements.size() && elements.get(i).getOperator() == CigarOperator.SOFT_CLIP ) {
            size += elements.get(i).getLength();
            i++;
        }

        return size;
    }

    private int calculateAlignmentStartShift(final Cigar oldCigar, final Cigar newCigar) {
        final int newShift = calcHardSoftOffset(newCigar);
        final int oldShift = calcHardSoftOffset(oldCigar);
        return newShift - oldShift;
    }

    private int calculateHardClippingAlignmentShift(final CigarElement cigarElement, final int clippedLength) {
        // Insertions should be discounted from the total hard clip count
        if (cigarElement.getOperator() == CigarOperator.INSERTION) {
            return -clippedLength;
        }// Deletions and Ns should be added to the total hard clip count (because we want to maintain the original alignment start)
        else if (cigarElement.getOperator() == CigarOperator.DELETION || cigarElement.getOperator() == CigarOperator.SKIPPED_REGION) {
            return cigarElement.getLength();
        }

        // There is no shift if we are not clipping an indel
        return 0;
    }

    private static final class CigarShift {
        private final Cigar cigar;
        private final int shiftFromStart;
        private final int shiftFromEnd;

        private CigarShift(final Cigar cigar, final int shiftFromStart, final int shiftFromEnd) {
            this.cigar = cigar;
            this.shiftFromStart = shiftFromStart;
            this.shiftFromEnd = shiftFromEnd;
        }
    }
}