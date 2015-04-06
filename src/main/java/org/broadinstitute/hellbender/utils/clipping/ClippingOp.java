package org.broadinstitute.hellbender.utils.clipping;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.Iterator;
import java.util.List;
import java.util.Stack;
import java.util.Vector;

import static org.broadinstitute.hellbender.utils.read.ReadUtils.*;

/**
 * Represents a clip on a read.  It has a type (see the enum) along with a start and stop in the bases
 * of the read, plus an option extraInfo (useful for carrying info where needed).
 * <p/>
 * Also holds the critical apply function that actually execute the clipping operation on a provided read,
 * according to the wishes of the supplied ClippingAlgorithm enum.
 */
public class ClippingOp {
    public final int start, stop; // inclusive

    public ClippingOp(int start, int stop) {
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
    public SAMRecord apply(ClippingRepresentation algorithm, SAMRecord originalRead) {
        SAMRecord read = ReadUtils.clone(originalRead);
        byte[] quals = read.getBaseQualities();
        byte[] bases = read.getReadBases();
        byte[] newBases = new byte[bases.length];
        byte[] newQuals = new byte[quals.length];

        switch (algorithm) {
            // important note:
            //   it's not safe to call read.getReadBases()[i] = 'N' or read.getBaseQualities()[i] = 0
            //   because you're not guaranteed to get a pointer to the actual array of bytes in the GATKSAMRecord
            case WRITE_NS:
                for (int i = 0; i < bases.length; i++) {
                    if (i >= start && i <= stop) {
                        newBases[i] = 'N';
                    }
                    else {
                        newBases[i] = bases[i];
                    }
                }
                read.setReadBases(newBases);
                break;
            case WRITE_Q0S:
                for (int i = 0; i < quals.length; i++) {
                    if (i >= start && i <= stop) {
                        newQuals[i] = 0;
                    }
                    else {
                        newQuals[i] = quals[i];
                    }
                }
                read.setBaseQualities(newQuals);
                break;
            case WRITE_NS_Q0S:
                for (int i = 0; i < bases.length; i++) {
                    if (i >= start && i <= stop) {
                        newQuals[i] = 0;
                        newBases[i] = 'N';
                    }
                    else {
                        newQuals[i] = quals[i];
                        newBases[i] = bases[i];
                    }
                }
                read.setBaseQualities(newBases);
                read.setReadBases(newBases);
                break;
            case HARDCLIP_BASES:
                read = hardClip(read, start, stop);
                break;

            case SOFTCLIP_BASES:
                if (read.getReadUnmappedFlag()) {
                    // we can't process unmapped reads
                    throw new UserException("Read Clipper cannot soft clip unmapped reads");
                }

                //System.out.printf("%d %d %d%n", stop, start, read.getReadLength());
                int myStop = stop;
                if ((stop + 1 - start) == read.getReadLength()) {
                    // BAM representation issue -- we can't SOFTCLIP away all bases in a read, just leave it alone
                    //Walker.logger.info(String.format("Warning, read %s has all bases clip but this can't be represented with SOFTCLIP_BASES, just leaving it alone", read.getReadName()));
                    //break;
                    myStop--; // just decrement stop
                }

                if (start > 0 && myStop != read.getReadLength() - 1)
                    throw new RuntimeException(String.format("Cannot apply soft clipping operator to the middle of a read: %s to be clipped at %d-%d", read.getReadName(), start, myStop));

                Cigar oldCigar = read.getCigar();

                int scLeft = 0, scRight = read.getReadLength();
                if (start == 0)
                    scLeft = myStop + 1;
                else
                    scRight = start;

                Cigar newCigar = softClip(oldCigar, scLeft, scRight);
                read.setCigar(newCigar);

                int newClippedStart = getNewAlignmentStartOffset(newCigar, oldCigar);
                int newStart = read.getAlignmentStart() + newClippedStart;
                read.setAlignmentStart(newStart);

                break;

            case REVERT_SOFTCLIPPED_BASES:
                read = revertSoftClippedBases(read);
                break;

            default:
                throw new IllegalStateException("Unexpected Clipping operator type " + algorithm);
        }

        return read;
    }

    private SAMRecord revertSoftClippedBases(SAMRecord read) {
        SAMRecord unclipped = ReadUtils.clone(read);

        Cigar unclippedCigar = new Cigar();
        int matchesCount = 0;
        for (CigarElement element : read.getCigar().getCigarElements()) {
            if (element.getOperator() == CigarOperator.SOFT_CLIP || element.getOperator() == CigarOperator.MATCH_OR_MISMATCH)
                matchesCount += element.getLength();
            else if (matchesCount > 0) {
                unclippedCigar.add(new CigarElement(matchesCount, CigarOperator.MATCH_OR_MISMATCH));
                matchesCount = 0;
                unclippedCigar.add(element);
            } else
                unclippedCigar.add(element);
        }
        if (matchesCount > 0)
            unclippedCigar.add(new CigarElement(matchesCount, CigarOperator.MATCH_OR_MISMATCH));

        unclipped.setCigar(unclippedCigar);
        final int newStart = read.getAlignmentStart() + calculateAlignmentStartShift(read.getCigar(), unclippedCigar);
        unclipped.setAlignmentStart(newStart);

        if ( newStart <= 0 ) {
            // if the start of the unclipped read occurs before the contig,
            // we must hard clip away the bases since we cannot represent reads with
            // negative or 0 alignment start values in the SAMRecord (e.g., 0 means unaligned)
            return hardClip(unclipped, 0, - newStart);
        } else {
            return unclipped;
        }
    }

    /**
     * Given a cigar string, get the number of bases hard or soft clipped at the start
     */
    private int getNewAlignmentStartOffset(final Cigar __cigar, final Cigar __oldCigar) {
        int num = 0;
        for (CigarElement e : __cigar.getCigarElements()) {
            if (!e.getOperator().consumesReferenceBases()) {
                if (e.getOperator().consumesReadBases()) {
                    num += e.getLength();
                }
            } else {
                break;
            }
        }

        int oldNum = 0;
        int curReadCounter = 0;

        for (CigarElement e : __oldCigar.getCigarElements()) {
            int curRefLength = e.getLength();
            int curReadLength = e.getLength();
            if (!e.getOperator().consumesReadBases()) {
                curReadLength = 0;
            }

            boolean truncated = false;
            if (curReadCounter + curReadLength > num) {
                curReadLength = num - curReadCounter;
                curRefLength = num - curReadCounter;
                truncated = true;
            }

            if (!e.getOperator().consumesReferenceBases()) {
                curRefLength = 0;
            }

            curReadCounter += curReadLength;
            oldNum += curRefLength;

            if (curReadCounter > num || truncated) {
                break;
            }
        }

        return oldNum;
    }

    /**
     * Given a cigar string, soft clip up to startClipEnd and soft clip starting at endClipBegin
     */
    private Cigar softClip(final Cigar __cigar, final int __startClipEnd, final int __endClipBegin) {
        if (__endClipBegin <= __startClipEnd) {
            //whole thing should be soft clipped
            int cigarLength = 0;
            for (CigarElement e : __cigar.getCigarElements()) {
                cigarLength += e.getLength();
            }

            Cigar newCigar = new Cigar();
            newCigar.add(new CigarElement(cigarLength, CigarOperator.SOFT_CLIP));
            assert newCigar.isValid(null, -1) == null;
            return newCigar;
        }

        int curLength = 0;
        Vector<CigarElement> newElements = new Vector<>();
        for (CigarElement curElem : __cigar.getCigarElements()) {
            if (!curElem.getOperator().consumesReadBases()) {
                if (curElem.getOperator() == CigarOperator.HARD_CLIP || curLength > __startClipEnd && curLength < __endClipBegin) {
                    newElements.add(new CigarElement(curElem.getLength(), curElem.getOperator()));
                }
                continue;
            }

            int s = curLength;
            int e = curLength + curElem.getLength();
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
                assert midLength >= 0;
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

        Vector<CigarElement> finalNewElements = new Vector<CigarElement>();
        CigarElement lastElement = null;
        for (CigarElement elem : newElements) {
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

        Cigar newCigar = new Cigar(finalNewElements);
        assert newCigar.isValid(null, -1) == null;
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
     * @param read a non-null read
     * @param start a start >= 0 and < read.length
     * @param stop a stop >= 0 and < read.length.
     * @return a cloned version of read that has been properly trimmed down
     */
    private SAMRecord hardClip(SAMRecord read, int start, int stop) {

        // If the read is unmapped there is no Cigar string and neither should we create a new cigar string
        final CigarShift cigarShift = (read.getReadUnmappedFlag()) ? new CigarShift(new Cigar(), 0, 0) : hardClipCigar(read.getCigar(), start, stop);

        // the cigar may force a shift left or right (or both) in case we are left with insertions
        // starting or ending the read after applying the hard clip on start/stop.
        final int newLength = read.getReadLength() - (stop - start + 1) - cigarShift.shiftFromStart - cigarShift.shiftFromEnd;
        final byte[] newBases = new byte[newLength];
        final byte[] newQuals = new byte[newLength];
        final int copyStart = (start == 0) ? stop + 1 + cigarShift.shiftFromStart : cigarShift.shiftFromStart;

        System.arraycopy(read.getReadBases(), copyStart, newBases, 0, newLength);
        System.arraycopy(read.getBaseQualities(), copyStart, newQuals, 0, newLength);

        final SAMRecord hardClippedRead = ReadUtils.clone(read);

        hardClippedRead.setBaseQualities(newQuals);
        hardClippedRead.setReadBases(newBases);
        hardClippedRead.setCigar(cigarShift.cigar);
        if (start == 0)
            hardClippedRead.setAlignmentStart(read.getAlignmentStart() + calculateAlignmentStartShift(read.getCigar(), cigarShift.cigar));

        if (hasBaseIndelQualities(read)) {
            final byte[] newBaseInsertionQuals = new byte[newLength];
            final byte[] newBaseDeletionQuals = new byte[newLength];
            System.arraycopy(getBaseInsertionQualities(read), copyStart, newBaseInsertionQuals, 0, newLength);
            System.arraycopy(getBaseDeletionQualities(read), copyStart, newBaseDeletionQuals, 0, newLength);
            setInsertionBaseQualities(hardClippedRead, newBaseInsertionQuals);
            setDeletionBaseQualities(hardClippedRead, newBaseDeletionQuals);
        }

        return hardClippedRead;

    }

    private CigarShift hardClipCigar(Cigar cigar, int start, int stop) {
        Cigar newCigar = new Cigar();
        int index = 0;
        int totalHardClipCount = stop - start + 1;
        int alignmentShift = 0; // caused by hard clipping deletions

        // hard clip the beginning of the cigar string
        if (start == 0) {
            Iterator<CigarElement> cigarElementIterator = cigar.getCigarElements().iterator();
            CigarElement cigarElement = cigarElementIterator.next();
            // Skip all leading hard clips
            while (cigarElement.getOperator() == CigarOperator.HARD_CLIP) {
                totalHardClipCount += cigarElement.getLength();
                if (cigarElementIterator.hasNext())
                    cigarElement = cigarElementIterator.next();
                else
                    throw new GATKException("Read is entirely hardclipped, shouldn't be trying to clip it's cigar string");
            }
            // keep clipping until we hit stop
            while (index <= stop) {
                int shift = 0;
                if (cigarElement.getOperator().consumesReadBases())
                    shift = cigarElement.getLength();

                // we're still clipping or just finished perfectly
                if (index + shift == stop + 1) {
                    alignmentShift += calculateHardClippingAlignmentShift(cigarElement, cigarElement.getLength());
                    newCigar.add(new CigarElement(totalHardClipCount + alignmentShift, CigarOperator.HARD_CLIP));
                }
                // element goes beyond what we need to clip
                else if (index + shift > stop + 1) {
                    int elementLengthAfterChopping = cigarElement.getLength() - (stop - index + 1);
                    alignmentShift += calculateHardClippingAlignmentShift(cigarElement, stop - index + 1);
                    newCigar.add(new CigarElement(totalHardClipCount + alignmentShift, CigarOperator.HARD_CLIP));
                    newCigar.add(new CigarElement(elementLengthAfterChopping, cigarElement.getOperator()));
                }
                index += shift;
                alignmentShift += calculateHardClippingAlignmentShift(cigarElement, shift);

                if (index <= stop && cigarElementIterator.hasNext())
                    cigarElement = cigarElementIterator.next();
                else
                    break;
            }

            // add the remaining cigar elements
            while (cigarElementIterator.hasNext()) {
                cigarElement = cigarElementIterator.next();
                newCigar.add(new CigarElement(cigarElement.getLength(), cigarElement.getOperator()));
            }
        }

        // hard clip the end of the cigar string
        else {
            Iterator<CigarElement> cigarElementIterator = cigar.getCigarElements().iterator();
            CigarElement cigarElement = cigarElementIterator.next();

            // Keep marching on until we find the start
            while (index < start) {
                int shift = 0;
                if (cigarElement.getOperator().consumesReadBases())
                    shift = cigarElement.getLength();

                // we haven't gotten to the start yet, keep everything as is.
                if (index + shift < start)
                    newCigar.add(new CigarElement(cigarElement.getLength(), cigarElement.getOperator()));

                    // element goes beyond our clip starting position
                else {
                    int elementLengthAfterChopping = start - index;
                    alignmentShift += calculateHardClippingAlignmentShift(cigarElement, cigarElement.getLength() - (start - index));

                    // if this last element is a HARD CLIP operator, just merge it with our hard clip operator to be added later
                    if (cigarElement.getOperator() == CigarOperator.HARD_CLIP)
                        totalHardClipCount += elementLengthAfterChopping;
                        // otherwise, maintain what's left of this last operator
                    else
                        newCigar.add(new CigarElement(elementLengthAfterChopping, cigarElement.getOperator()));
                }
                index += shift;
                if (index < start && cigarElementIterator.hasNext())
                    cigarElement = cigarElementIterator.next();
                else
                    break;
            }

            // check if we are hard clipping indels
            while (cigarElementIterator.hasNext()) {
                cigarElement = cigarElementIterator.next();
                alignmentShift += calculateHardClippingAlignmentShift(cigarElement, cigarElement.getLength());

                // if the read had a HardClip operator in the end, combine it with the Hard Clip we are adding
                if (cigarElement.getOperator() == CigarOperator.HARD_CLIP)
                    totalHardClipCount += cigarElement.getLength();
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

        for (final CigarElement cigarElement : cigar.getCigarElements())
            cigarStack.push(cigarElement);

        for (int i = 1; i <= 2; i++) {
            int shift = 0;
            int totalHardClip = 0;
            boolean readHasStarted = false;
            boolean addedHardClips = false;

            while (!cigarStack.empty()) {
                CigarElement cigarElement = cigarStack.pop();

                if (!readHasStarted &&
                        cigarElement.getOperator() != CigarOperator.DELETION &&
                        cigarElement.getOperator() != CigarOperator.SKIPPED_REGION &&
                        cigarElement.getOperator() != CigarOperator.HARD_CLIP)
                    readHasStarted = true;

                else if (!readHasStarted && cigarElement.getOperator() == CigarOperator.HARD_CLIP)
                    totalHardClip += cigarElement.getLength();

                else if (!readHasStarted && cigarElement.getOperator() == CigarOperator.DELETION)
                    totalHardClip += cigarElement.getLength();

                else if (!readHasStarted && cigarElement.getOperator() == CigarOperator.SKIPPED_REGION)
                    totalHardClip += cigarElement.getLength();

                if (readHasStarted) {
                    if (i == 1) {
                        if (!addedHardClips) {
                            if (totalHardClip > 0)
                                inverseCigarStack.push(new CigarElement(totalHardClip, CigarOperator.HARD_CLIP));
                            addedHardClips = true;
                        }
                        inverseCigarStack.push(cigarElement);
                    } else {
                        if (!addedHardClips) {
                            if (totalHardClip > 0)
                                cleanCigar.add(new CigarElement(totalHardClip, CigarOperator.HARD_CLIP));
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

    private int calculateAlignmentStartShift(Cigar oldCigar, Cigar newCigar) {
        final int newShift = calcHardSoftOffset(newCigar);
        final int oldShift = calcHardSoftOffset(oldCigar);
        return newShift - oldShift;
    }

    private int calculateHardClippingAlignmentShift(CigarElement cigarElement, int clippedLength) {
        // Insertions should be discounted from the total hard clip count
        if (cigarElement.getOperator() == CigarOperator.INSERTION)
            return -clippedLength;

            // Deletions and Ns should be added to the total hard clip count (because we want to maintain the original alignment start)
        else if (cigarElement.getOperator() == CigarOperator.DELETION || cigarElement.getOperator() == CigarOperator.SKIPPED_REGION)
            return cigarElement.getLength();

        // There is no shift if we are not clipping an indel
        return 0;
    }

    private static class CigarShift {
        private Cigar cigar;
        private int shiftFromStart;
        private int shiftFromEnd;

        private CigarShift(Cigar cigar, int shiftFromStart, int shiftFromEnd) {
            this.cigar = cigar;
            this.shiftFromStart = shiftFromStart;
            this.shiftFromEnd = shiftFromEnd;
        }
    }
}