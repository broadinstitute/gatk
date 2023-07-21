package org.broadinstitute.hellbender.utils.clipping;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

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
     *  @param algorithm    clipping algorithm to use
     * @param originalRead the read to be clipped
     */
    GATKRead apply(final ClippingRepresentation algorithm, final GATKRead originalRead) {
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
                return applySoftClipBases(originalRead.copy());
            }

            case REVERT_SOFTCLIPPED_BASES: {
                return applyRevertSoftClippedBases(originalRead.copy());
            }

            default: {
                throw new IllegalStateException("Unexpected Clipping operator type " + algorithm);
            }
        }
    }

    private GATKRead applySoftClipBases(final GATKRead readCopied) {
        Utils.validateArg(!readCopied.isUnmapped(), "Read Clipper cannot soft clip unmapped reads");

        if (readCopied.getLength() <= 2) {  // see GATK issue #2022 -- we can't soft-clip all bases
            return readCopied;
        }

        // see GATK issue #2022 -- we can't soft-clip all bases
        final int myStop = Math.min(stop, start + readCopied.getLength() - 2);
        Utils.validate(start <= 0 || myStop == readCopied.getLength() - 1, () -> String.format("Cannot apply soft clipping operator to the middle of a read: %s to be clipped at %d-%d", readCopied.getName(), start, myStop));

        final Cigar oldCigar = readCopied.getCigar();
        final Cigar newCigar = CigarUtils.clipCigar(oldCigar, start, myStop + 1, CigarOperator.SOFT_CLIP);
        readCopied.setCigar(newCigar);

        final int alignmentStartShift = start == 0 ? CigarUtils.alignmentStartShift(oldCigar, stop + 1) : 0;
        final int newStart = readCopied.getStart() + alignmentStartShift;
        readCopied.setPosition(readCopied.getContig(), newStart);
        return readCopied;
    }

    private void applyWriteQ0s(final GATKRead readCopied) {
        final byte[] newQuals = readCopied.getBaseQualities(); //this makes a copy so we can modify in place
        overwriteFromStartToStop(newQuals, (byte) 0);
        readCopied.setBaseQualities(newQuals);
    }

    private void applyWriteNs(final GATKRead readCopied) {
        final byte[] newBases = readCopied.getBases();       //this makes a copy so we can modify in place
        overwriteFromStartToStop(newBases, Nucleotide.N.encodeAsByte());
        readCopied.setBases(newBases);
    }

    private void overwriteFromStartToStop(final byte[] arr, final byte newVal) {
        Arrays.fill(arr, start, Math.min(arr.length, stop + 1), newVal);
    }

    private GATKRead applyRevertSoftClippedBases(final GATKRead read) {
        final Cigar originalCigar = read.getCigar();
        final List<CigarElement> originalElements = originalCigar.getCigarElements();
        if (originalElements.isEmpty() || !(originalElements.get(0).getOperator().isClipping() || originalElements.get(originalElements.size() - 1).getOperator().isClipping())) {
            return read;
        }
        GATKRead unclipped = read.copy();
        final Cigar unclippedCigar = CigarUtils.revertSoftClips(originalCigar);
        unclipped.setCigar(unclippedCigar);

        final int newStart = read.getSoftStart();

        if (newStart <= 0) {
            // if the start of the unclipped read occurs before the contig,
            // we must hard clip away the bases since we cannot represent reads with
            // negative or 0 alignment start values in the SAMRecord (e.g., 0 means unaligned)

            // We cannot set the read to temporarily have a negative start position, as our Read
            // interface will not allow it. Instead, since we know that the final start position will
            // be 1 after the hard clip operation, set it to 1 explicitly. We have to set it twice:
            // once before the hard clip (to reset the alignment stop / read length in read implementations
            // that cache these values, such as SAMRecord), and again after the hard clip.
            unclipped.setPosition(unclipped.getContig(), 1);
            unclipped = applyHardClipBases(unclipped, 0, -newStart);

            // Reset the position to 1 again only if we didn't end up with an empty, unmapped read after hard clipping.
            // See https://github.com/broadinstitute/gatk/issues/3845
            if (!unclipped.isUnmapped()) {
                unclipped.setPosition(unclipped.getContig(), 1);
            }

            return unclipped;
        } else {
            unclipped.setPosition(unclipped.getContig(), newStart);
            return unclipped;
        }
    }

    /**
     * Hard clip bases from read, from start to stop in base coordinates
     * <p>
     * If start == 0, then we will clip from the front of the read, otherwise we clip
     * from the right.  If start == 0 and stop == 10, this would clip out the first
     * 10 bases of the read.
     * <p>
     * Note that this function works with reads with negative alignment starts, in order to
     * allow us to hardClip reads that have had their soft clips reverted and so might have
     * negative alignment starts
     * <p>
     * Works properly with reduced reads and insertion/deletion base qualities
     * <p>
     * Note: this method does not assume that the read is directly modifiable
     * and makes a copy of it.
     *
     * @param read  a non-null read
     * @param start a start >= 0 and < read.length
     * @param stop  a stop >= 0 and < read.length.
     * @return a cloned version of read that has been properly trimmed down (Could be an empty, unmapped read)
     */
    private GATKRead applyHardClipBases(final GATKRead read, final int start, final int stop) {
        final int newLength = read.getLength() - (stop - start + 1);

        // If the new read is going to be empty, return an empty read now. This avoids initializing the new
        // read with invalid values below in certain cases (such as a negative alignment start).
        // See https://github.com/broadinstitute/gatk/issues/3466
        if (newLength == 0) {
            return ReadUtils.emptyRead(read);
        }

        // If the read is unmapped there is no Cigar string and neither should we create a new cigar string

        final Cigar cigar = read.getCigar();//Get the cigar once to avoid multiple calls because each makes a copy of the cigar
        final Cigar newCigar = read.isUnmapped() ? new Cigar() : CigarUtils.clipCigar(cigar, start, stop + 1, CigarOperator.HARD_CLIP);

        final byte[] newBases = new byte[newLength];
        final byte[] newQuals = new byte[newLength];
        final int copyStart = (start == 0) ? stop + 1 : 0;

        System.arraycopy(read.getBases(), copyStart, newBases, 0, newLength);
        System.arraycopy(read.getBaseQualities(), copyStart, newQuals, 0, newLength);

        final GATKRead hardClippedRead = read.copy();
        hardClippedRead.hardClipAttributes(copyStart, newLength, read.getBasesNoCopy().length);

        hardClippedRead.setBaseQualities(newQuals);
        hardClippedRead.setBases(newBases);
        hardClippedRead.setCigar(newCigar);
        if (start == 0 && !read.isUnmapped()) {
            hardClippedRead.setPosition(read.getContig(), read.getStart() + CigarUtils.alignmentStartShift(cigar, stop + 1));
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
}