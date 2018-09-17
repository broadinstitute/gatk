package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;

import java.util.*;
import java.util.function.Function;


public final class AlignmentUtils {
    private static final EnumSet<CigarOperator> ALIGNED_TO_GENOME_OPERATORS = EnumSet.of(CigarOperator.M, CigarOperator.EQ, CigarOperator.X);
    private static final EnumSet<CigarOperator> ALIGNED_TO_GENOME_PLUS_SOFTCLIPS = EnumSet.of(CigarOperator.M, CigarOperator.EQ, CigarOperator.X, CigarOperator.S);
    public final static String HAPLOTYPE_TAG = "HC";
    public final static byte GAP_CHARACTER = (byte)'-';

    // cannot be instantiated
    private AlignmentUtils() { }

    /**
     * Does cigar start or end with a deletion operation?
     *
     * @param cigar a non-null cigar to test
     * @return true if the first or last operator of cigar is a D
     */
    public static boolean startsOrEndsWithInsertionOrDeletion(final Cigar cigar) {
        Utils.nonNull(cigar);

        if ( cigar.isEmpty() )
            return false;

        final CigarOperator first = cigar.getCigarElement(0).getOperator();
        final CigarOperator last = cigar.getCigarElement(cigar.numCigarElements()-1).getOperator();
        return first == CigarOperator.D || first == CigarOperator.I || last == CigarOperator.D || last == CigarOperator.I;
    }

    /**
     * Aligns reads the haplotype, and then projects this alignment of read -> hap onto the reference
     * via the alignment of haplotype (via its getCigar) method.
     *
     * @param originalRead the read we want to write aligned to the reference genome
     * @param haplotype the haplotype that the read should be aligned to, before aligning to the reference
     * @param referenceStart the start of the reference that haplotype is aligned to.  Provides global coordinate frame.
     * @param isInformative true if the read is differentially informative for one of the haplotypes
     *
     * @param aligner
     * @throws IllegalArgumentException if {@code originalRead} is {@code null} or {@code haplotype} is {@code null} or it
     *   does not have a Cigar or the {@code referenceStart} is invalid (less than 1).
     *
     * @return a GATKRead aligned to reference. Never {@code null}.
     */
    public static GATKRead createReadAlignedToRef(final GATKRead originalRead,
                                                  final Haplotype haplotype,
                                                  final Haplotype refHaplotype,
                                                  final int referenceStart,
                                                  final boolean isInformative,
                                                  final SmithWatermanAligner aligner) {
        Utils.nonNull(originalRead);
        Utils.nonNull(haplotype);
        Utils.nonNull(refHaplotype);
        Utils.nonNull(haplotype.getCigar());
        Utils.nonNull(aligner);
        if ( referenceStart < 1 ) { throw new IllegalArgumentException("reference start much be >= 1 but got " + referenceStart); }

        // compute the smith-waterman alignment of read -> haplotype
        final SmithWatermanAlignment swPairwiseAlignment = aligner.align(haplotype.getBases(), originalRead.getBases(), CigarUtils.ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS, SWOverhangStrategy.SOFTCLIP);
        if ( swPairwiseAlignment.getAlignmentOffset() == -1 ) {
            // sw can fail (reasons not clear) so if it happens just don't realign the read
            return originalRead;
        }

        final Cigar swCigar = consolidateCigar(swPairwiseAlignment.getCigar());

        // since we're modifying the read we need to clone it
        final GATKRead read = originalRead.copy();

        // only informative reads are given the haplotype tag to enhance visualization
        if ( isInformative ) {
            read.setAttribute(HAPLOTYPE_TAG, haplotype.hashCode());
        }

        // compute here the read starts w.r.t. the reference from the SW result and the hap -> ref cigar
        final Cigar extendedHaplotypeCigar = haplotype.getConsolidatedPaddedCigar(1000);
        final int readStartOnHaplotype = calcFirstBaseMatchingReferenceInCigar(extendedHaplotypeCigar, swPairwiseAlignment.getAlignmentOffset());
        final int readStartOnReference = referenceStart + haplotype.getAlignmentStartHapwrtRef() + readStartOnHaplotype;

        // compute the read -> ref alignment by mapping read -> hap -> ref from the
        // SW of read -> hap mapped through the given by hap -> ref
        final Cigar haplotypeToRef = trimCigarByBases(extendedHaplotypeCigar, swPairwiseAlignment.getAlignmentOffset(), extendedHaplotypeCigar.getReadLength() - 1);
        final Cigar readToRefCigarRaw = applyCigarToCigar(swCigar, haplotypeToRef);
        final Cigar readToRefCigarClean = cleanUpCigar(readToRefCigarRaw);
        final Cigar readToRefCigar = leftAlignIndel(readToRefCigarClean, refHaplotype.getBases(),
                originalRead.getBases(), readStartOnHaplotype, 0, true);

        final int leadingDeletions = readToRefCigarClean.getReferenceLength() - readToRefCigar.getReferenceLength();
        read.setPosition(read.getContig(), readStartOnReference + leadingDeletions);

        // the SW Cigar does not contain the hard clips of the original read
        final Cigar originalCigar = originalRead.getCigar();
        final CigarElement firstElement = originalCigar.getFirstCigarElement();
        final CigarElement lastElement = originalCigar.getLastCigarElement();
        final List<CigarElement> readToRefCigarElementsWithHardClips = new ArrayList<>();
        if (firstElement.getOperator() == CigarOperator.HARD_CLIP) {
            readToRefCigarElementsWithHardClips.add(firstElement);
        }
        readToRefCigarElementsWithHardClips.addAll(readToRefCigar.getCigarElements());
        if (lastElement.getOperator() == CigarOperator.HARD_CLIP) {
            readToRefCigarElementsWithHardClips.add(lastElement);
        }

        read.setCigar(new Cigar(readToRefCigarElementsWithHardClips));

        if ( readToRefCigar.getReadLength() != read.getLength() ) {
            throw new GATKException("Cigar " + readToRefCigar + " with read length " + readToRefCigar.getReadLength()
                    + " != read length " + read.getLength() + " for read " + read.toString() + "\nhapToRef " + haplotypeToRef + " length " + haplotypeToRef.getReadLength() + "/" + haplotypeToRef.getReferenceLength()
                    + "\nreadToHap " + swCigar + " length " + swCigar.getReadLength() + "/" + swCigar.getReferenceLength());
        }

        return read;
    }


    /**
     * Get the byte[] from bases that cover the reference interval refStart -> refEnd given the
     * alignment of bases to the reference (basesToRefCigar) and the start offset of the bases on the reference
     *
     * refStart and refEnd are 0 based offsets that we want to obtain.  In the client code, if the reference
     * bases start at position X and you want Y -> Z, refStart should be Y - X and refEnd should be Z - X.
     *
     * If refStart or refEnd would start or end the new bases within a deletion, this function will return null
     *
     * @param bases
     * @param refStart
     * @param refEnd
     * @param basesStartOnRef where does the bases array start w.r.t. the reference start?  For example, bases[0] of
     *                        could be at refStart == 0 if basesStartOnRef == 0, but it could just as easily be at
     *                        10 (meaning bases doesn't fully span the reference), which would be indicated by basesStartOnRef == 10.
     *                        It's not trivial to eliminate this parameter because it's tied up with the cigar
     * @param basesToRefCigar the cigar that maps the bases to the reference genome
     * @return a byte[] containing the bases covering this interval, or null if we would start or end within a deletion
     */
    public static byte[] getBasesCoveringRefInterval(final int refStart, final int refEnd, final byte[] bases, final int basesStartOnRef, final Cigar basesToRefCigar) {
        if ( refStart < 0 || refEnd < refStart ) throw new IllegalArgumentException("Bad start " + refStart + " and/or stop " + refEnd);
        if ( basesStartOnRef < 0 ) throw new IllegalArgumentException("BasesStartOnRef must be >= 0 but got " + basesStartOnRef);
        Utils.nonNull( bases );
        Utils.nonNull( basesToRefCigar );
        if ( bases.length != basesToRefCigar.getReadLength() ) throw new IllegalArgumentException("Mismatch in length between reference bases " + bases.length + " and cigar length " + basesToRefCigar);

        int refPos = basesStartOnRef;
        int basesPos = 0;
        int basesStart = -1;
        int basesStop = -1;
        boolean done = false;

        for ( int iii = 0; ! done && iii < basesToRefCigar.numCigarElements(); iii++ ) {
            final CigarElement ce = basesToRefCigar.getCigarElement(iii);
            switch ( ce.getOperator() ) {
                case I:
                    basesPos += ce.getLength();
                    break;
                case M: case X: case EQ:
                    for ( int i = 0; i < ce.getLength(); i++ ) {
                        if ( refPos == refStart )
                            basesStart = basesPos;
                        if ( refPos == refEnd ) {
                            basesStop = basesPos;
                            done = true;
                            break;
                        }
                        refPos++;
                        basesPos++;
                    }
                    break;
                case D:
                    for ( int i = 0; i < ce.getLength(); i++ ) {
                        if ( refPos == refEnd || refPos == refStart ) {
                            // if we ever reach a ref position that is either a start or an end, we fail
                            return null;
                        }
                        refPos++;
                    }
                    break;
                default:
                    throw new IllegalStateException("Unsupported operator " + ce);
            }
        }

        if ( basesStart == -1 || basesStop == -1 )
            throw new IllegalStateException("Never found start " + basesStart + " or stop " + basesStop + " given cigar " + basesToRefCigar);

        return Arrays.copyOfRange(bases, basesStart, basesStop + 1);
    }

    public static byte[] getBasesAlignedOneToOne(final GATKRead read) {
        return getSequenceAlignedOneToOne(read, r -> r.getBasesNoCopy(), GAP_CHARACTER);
    }

    public static byte[] getBaseQualsAlignedOneToOne(final GATKRead read) {
        return getSequenceAlignedOneToOne(read, r -> r.getBaseQualitiesNoCopy(), (byte)0);
    }

    public static byte[] getSequenceAlignedOneToOne(final GATKRead read, final Function<GATKRead, byte[]> bytesProvider, final byte padWith) {
        Utils.nonNull(read);
        Utils.nonNull(bytesProvider);
        final Cigar cigar = read.getCigar();
        final byte[] sequence = bytesProvider.apply(read);

        if (!cigar.containsOperator(CigarOperator.DELETION) && !cigar.containsOperator(CigarOperator.INSERTION)) {
            return sequence;
        }
        else {
            final byte[] paddedBases = new byte[CigarUtils.countRefBasesIncludingSoftClips(read, 0, cigar.numCigarElements())];
            int literalPos = 0;
            int paddedPos = 0;
            for ( int i = 0; i < cigar.numCigarElements(); i++ ) {
                final CigarElement ce = cigar.getCigarElement(i);
                final CigarOperator co = ce.getOperator();
                if (co.consumesReadBases()) {
                    if (!co.consumesReferenceBases()) {
                        literalPos += ce.getLength();  //skip inserted bases
                    }
                    else {
                        System.arraycopy(sequence, literalPos, paddedBases, paddedPos, ce.getLength());
                        literalPos += ce.getLength();
                        paddedPos += ce.getLength();
                    }
                }
                else if (co.consumesReferenceBases()) {
                    for ( int j = 0; j < ce.getLength(); j++ ) {  //pad deleted bases
                        paddedBases[paddedPos] = padWith;
                        paddedPos++;
                    }
                }
            }
            return paddedBases;
        }
    }

    /**
     * Get the number of bases at which refSeq and readSeq differ, given their alignment
     *
     * @param cigar the alignment of readSeq to refSeq
     * @param refSeq the bases of the reference sequence
     * @param readSeq the bases of the read sequence
     * @return the number of bases that differ between refSeq and readSeq
     */
    @SuppressWarnings("fallthrough")
    public static int calcNumDifferentBases(final Cigar cigar, final byte[] refSeq, final byte[] readSeq) {
        int refIndex = 0, readIdx = 0, delta = 0;

        for (final CigarElement ce : cigar.getCigarElements()) {
            final int elementLength = ce.getLength();
            switch (ce.getOperator()) {
                case X:case EQ:case M:
                    for (int j = 0; j < elementLength; j++, refIndex++, readIdx++)
                        delta += refSeq[refIndex] != readSeq[readIdx] ? 1 : 0;
                    break;
                case I:
                    delta += elementLength;
                case S:
                    readIdx += elementLength;
                    break;
                case D:
                    delta += elementLength;
                case N:
                    refIndex += elementLength;
                    break;
                case H:
                case P:
                    break;
                default:
                    throw new GATKException("The " + ce.getOperator() + " cigar element is not currently supported");
            }
        }

        return delta;
    }

    public static class MismatchCount {
        public int numMismatches = 0;
        public long mismatchQualities = 0;
    }

    // todo -- this code and mismatchesInRefWindow should be combined and optimized into a single
    // todo -- high performance implementation.  We can do a lot better than this right now
    /**
     * @see #getMismatchCount(GATKRead, byte[], int, int, int) with startOnRead == 0 and nReadBases == read.getReadLength()
     */
    public static MismatchCount getMismatchCount(GATKRead r, byte[] refSeq, int refIndex) {
        return getMismatchCount(r, refSeq, refIndex, 0, r.getLength());
    }

    /**
     * Count how many bases mismatch the reference.  Indels are not considered mismatching.
     *
     * @param r                   the sam record to check against
     * @param refSeq              the byte array representing the reference sequence
     * @param refIndex            the index in the reference byte array of the read's first base (the reference index
     *                            is matching the alignment start, there may be tons of soft-clipped bases before/after
     *                            that so it's wrong to compare with getLength() here.).  Note that refIndex is
     *                            zero based, not 1 based
     * @param startOnRead         the index in the read's bases from which we start counting
     * @param nReadBases          the number of bases after (but including) startOnRead that we check
     * @return non-null object representing the mismatch count
     */
    @SuppressWarnings("fallthrough")
    public static MismatchCount getMismatchCount(GATKRead r, byte[] refSeq, int refIndex, int startOnRead, int nReadBases) {
        Utils.nonNull( r );
        Utils.nonNull( refSeq );
        if ( refIndex < 0 ) throw new IllegalArgumentException("attempting to calculate the mismatch count with a reference index that is negative");
        if ( startOnRead < 0 ) throw new IllegalArgumentException("attempting to calculate the mismatch count with a read start that is negative");
        if ( nReadBases < 0 ) throw new IllegalArgumentException("attempting to calculate the mismatch count for a negative number of read bases");
        if ( refSeq.length - refIndex < (r.getEnd() - r.getStart()) )
            throw new IllegalArgumentException("attempting to calculate the mismatch count against a reference string that is smaller than the read");

        MismatchCount mc = new MismatchCount();

        int readIdx = 0;
        final int endOnRead = startOnRead + nReadBases - 1; // index of the last base on read we want to count (note we are including soft-clipped bases with this math)
        final byte[] readSeq = r.getBases();
        final Cigar c = r.getCigar();
        final byte[] readQuals = r.getBaseQualities();
        for (final CigarElement ce : c.getCigarElements()) {

            if (readIdx > endOnRead)
                break;

            final int elementLength = ce.getLength();
            switch (ce.getOperator()) {
                case X:
                    mc.numMismatches += elementLength;
                    for (int j = 0; j < elementLength; j++)
                        mc.mismatchQualities += readQuals[readIdx+j];
                case EQ:
                    refIndex += elementLength;
                    readIdx += elementLength;
                    break;
                case M:
                    for (int j = 0; j < elementLength; j++, refIndex++, readIdx++) {
                        if (refIndex >= refSeq.length)
                            continue;                      // TODO : It should never happen, we should throw exception here
                        if (readIdx < startOnRead) continue;
                        if (readIdx > endOnRead) break;
                        byte refChr = refSeq[refIndex];
                        byte readChr = readSeq[readIdx];
                        // Note: we need to count X/N's as mismatches because that's what SAM requires
                        //if ( BaseUtils.simpleBaseToBaseIndex(readChr) == -1 ||
                        //     BaseUtils.simpleBaseToBaseIndex(refChr)  == -1 )
                        //    continue; // do not count Ns/Xs/etc ?
                        if (readChr != refChr) {
                            mc.numMismatches++;
                            mc.mismatchQualities += readQuals[readIdx];
                        }
                    }
                    break;
                case I:
                case S:
                    readIdx += elementLength;
                    break;
                case D:
                case N:
                    refIndex += elementLength;
                    break;
                case H:
                case P:
                    break;
                default:
                    throw new GATKException("The " + ce.getOperator() + " cigar element is not currently supported");
            }

        }
        return mc;
    }

    /**
     * Returns number of alignment blocks (continuous stretches of aligned bases) in the specified alignment.
     * This method follows closely the SAMRecord::getAlignmentBlocks() implemented in samtools library, but
     * it only counts blocks without actually allocating and filling the list of blocks themselves. Hence, this method is
     * a much more efficient alternative to r.getAlignmentBlocks.size() in the situations when this number is all that is needed.
     * Formally, this method simply returns the number of M elements in the cigar.
     *
     * @param r alignment
     * @return number of continuous alignment blocks (i.e. 'M' elements of the cigar; all indel and clipping elements are ignored).
     */
    public static int getNumAlignmentBlocks(final GATKRead r) {
        Utils.nonNull( r );
        final Cigar cigar = r.getCigar();
        if (cigar == null) return 0;

        int n = 0;
        for (final CigarElement e : cigar.getCigarElements()) {
            if (ALIGNED_TO_GENOME_OPERATORS.contains(e.getOperator()))
                n++;
        }

        return n;
    }


    /**
     * Get the number of bases aligned to the genome, including soft clips
     *
     * If read is not mapped (i.e., doesn't have a cigar) returns 0
     *
     * @param r a non-null Read
     * @return the number of bases aligned to the genome in R, including soft clipped bases
     */
    public static int getNumAlignedBasesCountingSoftClips(final GATKRead r) {
        int n = 0;
        final Cigar cigar = r.getCigar();
        if (cigar == null) return 0;

        for (final CigarElement e : cigar.getCigarElements())
            if (ALIGNED_TO_GENOME_PLUS_SOFTCLIPS.contains(e.getOperator()))
                n += e.getLength();

        return n;
    }

    /**
     * Count the number of bases hard clipped from read
     *
     * If read's cigar is null, return 0
     *
     * @param r a non-null read
     * @return a positive integer
     */
    public static int getNumHardClippedBases(final GATKRead r) {
        if ( r == null ) throw new IllegalArgumentException("Read cannot be null");

        int n = 0;
        final Cigar cigar = r.getCigar();
        if (cigar == null) return 0;

        for (final CigarElement e : cigar.getCigarElements())
            if (e.getOperator() == CigarOperator.H)
                n += e.getLength();

        return n;
    }

    /**
     * Calculate the number of bases that are soft clipped in read with quality score greater than threshold
     *
     * Handles the case where the cigar is null (i.e., the read is unmapped), returning 0
     *
     * @param read a non-null read.
     * @param qualThreshold consider bases with quals > this value as high quality.  Must be >= 0
     * @return positive integer
     */
    public static int calcNumHighQualitySoftClips( final GATKRead read, final byte qualThreshold ) {
        if ( read == null ) throw new IllegalArgumentException("Read cannot be null");
        if ( qualThreshold < 0 ) throw new IllegalArgumentException("Expected qualThreshold to be a positive byte but saw " + qualThreshold);

        if ( read.getCigar() == null ) // the read is unmapped
            return 0;

        final byte[] qual = read.getBaseQualities();

        int numHQSoftClips = 0;
        int alignPos = 0;
        for ( final CigarElement ce : read.getCigarElements() ) {
            final int elementLength = ce.getLength();

            switch( ce.getOperator() ) {
                case S:
                    for( int jjj = 0; jjj < elementLength; jjj++ ) {
                        if( qual[alignPos++] > qualThreshold ) { numHQSoftClips++; }
                    }
                    break;
                case M: case I: case EQ: case X:
                    alignPos += elementLength;
                    break;
                case H: case P: case D: case N:
                    break;
                default:
                    throw new IllegalStateException("Unsupported cigar operator: " + ce.getOperator());
            }
        }

        return numHQSoftClips;
    }

    public static int calcAlignmentByteArrayOffset(final Cigar cigar, final PileupElement pileupElement, final int alignmentStart, final int refLocus) {
        return calcAlignmentByteArrayOffset( cigar, pileupElement.getOffset(), pileupElement.isDeletion(), alignmentStart, refLocus );
    }

    /**
     * Calculate the index into the read's bases of the beginning of the encompassing cigar element for a given cigar and offset
     *
     * @param cigar            the read's CIGAR -- cannot be null
     * @param offset           the offset to use for the calculation or -1 if in the middle of a deletion
     * @param isDeletion       are we in the middle of a deletion?
     * @param alignmentStart   the alignment start of the read
     * @param refLocus         the reference position of the offset
     * @return a non-negative int index
     */
    public static int calcAlignmentByteArrayOffset(final Cigar cigar, final int offset, final boolean isDeletion, final int alignmentStart, final int refLocus) {
        if ( cigar == null ) throw new IllegalArgumentException("attempting to find the alignment position from a CIGAR that is null");
        if ( offset < -1 ) throw new IllegalArgumentException("attempting to find the alignment position with an offset that is negative (and not -1)");
        if ( alignmentStart < 0 ) throw new IllegalArgumentException("attempting to find the alignment position from an alignment start that is negative");
        if ( refLocus < 0 ) throw new IllegalArgumentException("attempting to find the alignment position from a reference position that is negative");
        if ( offset >= cigar.getReadLength() ) throw new IllegalArgumentException("attempting to find the alignment position of an offset than is larger than the read length");

        int pileupOffset = offset;

        // Reassign the offset if we are in the middle of a deletion because of the modified representation of the read bases
        if (isDeletion) {
            pileupOffset = refLocus - alignmentStart;
            final CigarElement ce = cigar.getCigarElement(0);
            if (ce.getOperator() == CigarOperator.S) {
                pileupOffset += ce.getLength();
            }
        }

        int pos = 0;
        int alignmentPos = 0;

        for (int iii = 0; iii < cigar.numCigarElements(); iii++) {
            final CigarElement ce = cigar.getCigarElement(iii);
            final int elementLength = ce.getLength();

            switch (ce.getOperator()) {
                case I:
                case S: // TODO -- I don't think that soft clips should be treated the same as inserted bases here. Investigation needed.
                    pos += elementLength;
                    if (pos >= pileupOffset) {
                        return alignmentPos;
                    }
                    break;
                case D:
                    if (!isDeletion) {
                        alignmentPos += elementLength;
                    } else {
                        if (pos + elementLength - 1 >= pileupOffset) {
                            return alignmentPos + (pileupOffset - pos);
                        } else {
                            pos += elementLength;
                            alignmentPos += elementLength;
                        }
                    }
                    break;
                case M:
                case EQ:
                case X:
                    if (pos + elementLength - 1 >= pileupOffset) {
                        return alignmentPos + (pileupOffset - pos);
                    } else {
                        pos += elementLength;
                        alignmentPos += elementLength;
                    }
                    break;
                case H:
                case P:
                case N:
                    break;
                default:
                    throw new GATKException("Unsupported cigar operator: " + ce.getOperator());
            }
        }

        return alignmentPos;
    }

    /**
     * Is the offset inside a deletion?
     *
     * @param cigar         the read's CIGAR -- cannot be null
     * @param offset        the offset into the CIGAR
     * @return true if the offset is inside a deletion, false otherwise
     */
    public static boolean isInsideDeletion(final Cigar cigar, final int offset) {
        Utils.nonNull(cigar);
        if ( offset < 0 ) return false;

        // pos counts read bases
        int pos = 0;
        int prevPos = 0;

        for (final CigarElement ce : cigar.getCigarElements()) {

            switch (ce.getOperator()) {
                case I:
                case S:
                case D:
                case M:
                case EQ:
                case X:
                    prevPos = pos;
                    pos += ce.getLength();
                    break;
                case H:
                case P:
                case N:
                    break;
                default:
                    throw new GATKException("Unsupported cigar operator: " + ce.getOperator());
            }

            // Is the offset inside a deletion?
            if ( prevPos < offset && pos >= offset && ce.getOperator() == CigarOperator.D ) {
                return true;

            }
        }

        return false;
    }

    /**
     * Generate an array of bases for just those that are aligned to the reference (i.e. no clips or insertions)
     *
     * @param cigar            the read's CIGAR -- cannot be null
     * @param read             the read's base array
     * @return a non-null array of bases (bytes)
     */
    @SuppressWarnings("fallthrough")
    public static byte[] readToAlignmentByteArray(final Cigar cigar, final byte[] read) {
        Utils.nonNull(cigar);
        Utils.nonNull(read);

        final int alignmentLength = cigar.getReferenceLength();
        final byte[] alignment = new byte[alignmentLength];
        int alignPos = 0;
        int readPos = 0;
        for (int iii = 0; iii < cigar.numCigarElements(); iii++) {

            final CigarElement ce = cigar.getCigarElement(iii);
            final int elementLength = ce.getLength();

            switch (ce.getOperator()) {
                case I:
                    if (alignPos > 0) {
                        final int prevPos = alignPos - 1;
                        if (alignment[prevPos] == BaseUtils.Base.A.base) {
                            alignment[prevPos] = PileupElement.A_FOLLOWED_BY_INSERTION_BASE;
                        } else if (alignment[prevPos] == BaseUtils.Base.C.base) {
                            alignment[prevPos] = PileupElement.C_FOLLOWED_BY_INSERTION_BASE;
                        } else if (alignment[prevPos] == BaseUtils.Base.T.base) {
                            alignment[prevPos] = PileupElement.T_FOLLOWED_BY_INSERTION_BASE;
                        } else if (alignment[prevPos] == BaseUtils.Base.G.base) {
                            alignment[prevPos] = PileupElement.G_FOLLOWED_BY_INSERTION_BASE;
                        }
                    }
                case S:
                    readPos += elementLength;
                    break;
                case D:
                case N:
                    for (int jjj = 0; jjj < elementLength; jjj++) {
                        alignment[alignPos++] = PileupElement.DELETION_BASE;
                    }
                    break;
                case M:
                case EQ:
                case X:
                    for (int jjj = 0; jjj < elementLength; jjj++) {
                        alignment[alignPos++] = read[readPos++];
                    }
                    break;
                case H:
                case P:
                    break;
                default:
                    throw new GATKException("Unsupported cigar operator: " + ce.getOperator());
            }
        }
        return alignment;
    }

    /**
     * Need a well-formed, consolidated Cigar string so that the left aligning code works properly.
     * For example, 1M1M1M1D2M1M --> 3M1D3M
     * If the given cigar is empty then the returned cigar will also be empty
     *
     * Note that this routine collapses cigar elements of size 0, so 2M0M => 2M
     *
     * @param c the cigar to consolidate
     * @return  a non-null cigar with consecutive matching operators merged into single operators.
     */
    public static Cigar consolidateCigar( final Cigar c ) {
        if ( c == null ) { throw new IllegalArgumentException("Cigar cannot be null"); }

        // fast check to determine if there's anything worth doing before we create new Cigar and actually do some work
        if ( ! needsConsolidation(c) )
            return c;

        final Cigar returnCigar = new Cigar();
        int sumLength = 0;
        CigarElement lastElement = null;

        for( final CigarElement cur : c.getCigarElements() ) {
            if ( cur.getLength() == 0 )
                continue; // don't add elements of 0 length

            if ( lastElement != null && lastElement.getOperator() != cur.getOperator() ) {
                returnCigar.add(new CigarElement(sumLength, lastElement.getOperator()));
                sumLength = 0;
            }

            sumLength += cur.getLength();
            lastElement = cur;
        }

        if ( sumLength > 0 ) {
            returnCigar.add(new CigarElement(sumLength, lastElement.getOperator()));
        }

        return returnCigar;
    }

    /**
     * Does the cigar C need to be consolidated?
     *
     * @param c a non-null cigar
     * @return true if so
     */
    private static boolean needsConsolidation(final Cigar c) {
        if ( c.numCigarElements() <= 1 )
            return false; // fast path for empty or single cigar

        CigarOperator lastOp = null;
        for( final CigarElement cur : c.getCigarElements() ) {
            if ( cur.getLength() == 0 || lastOp == cur.getOperator() )
                return true;
            lastOp = cur.getOperator();
        }

        return false;
    }

    /**
     * Takes the alignment of the read sequence <code>readSeq</code> to the reference sequence <code>refSeq</code>
     * starting at 0-based position <code>refIndex</code> on the <code>refSeq</code> and specified by its <code>cigar</code>.
     * The last argument <code>readIndex</code> specifies 0-based position on the read where the alignment described by the
     * <code>cigar</code> starts. Usually cigars specify alignments of the whole read to the ref, so that readIndex is normally 0.
     * Use non-zero readIndex only when the alignment cigar represents alignment of a part of the read. The refIndex in this case
     * should be the position where the alignment of that part of the read starts at. In other words, both refIndex and readIndex are
     * always the positions where the cigar starts on the ref and on the read, respectively.
     * <p/>
     * If the alignment has one or more indels, this method attempts to move them left across a stretch of repetitive bases.
     * For instance, if the original cigar specifies that (any) one AT is deleted from a repeat sequence TATATATA, the output
     * cigar will always mark the leftmost AT as deleted. If there is no indel in the original cigar or if the indel position
     * is determined unambiguously (i.e. inserted/deleted sequence is not repeated), the original cigar is returned.
     *
     * Note that currently we do not actually support the case where there is more than one indel in the alignment.  We will throw
     * an exception if there is -- unless the
     *
     * @param cigar     structure of the original alignment
     * @param refSeq    reference sequence the read is aligned to
     * @param readSeq   read sequence
     * @param refIndex  0-based alignment start position on ref
     * @param readIndex 0-based alignment start position on read
     * @param leftmostAllowedAlignment left align indel no further left than this index (0-based)
     * @param doNotThrowExceptionForMultipleIndels  if true we will not throw an exception if we encounter multiple indels in the alignment will instead will return the original cigar
     * @return a non-null cigar, in which the indels are guaranteed to be placed at the leftmost possible position across a repeat (if any)
     */

    public static Cigar leftAlignIndel(Cigar cigar, final byte[] refSeq, final byte[] readSeq, final int refIndex, final int readIndex, final int leftmostAllowedAlignment,final boolean doNotThrowExceptionForMultipleIndels) {
        ensureLeftAlignmentHasGoodArguments(cigar, refSeq, readSeq, refIndex, readIndex);

        final int numIndels = countIndelElements(cigar);
        if ( numIndels == 0 )
            return cigar;
        if ( numIndels == 1 )
            return leftAlignSingleIndel(cigar, refSeq, readSeq, refIndex, readIndex, leftmostAllowedAlignment,true);

        // if we got here then there is more than 1 indel in the alignment
        if ( doNotThrowExceptionForMultipleIndels )
            return cigar;

        throw new UnsupportedOperationException("attempting to left align a CIGAR that has more than 1 indel in its alignment but this functionality has not been implemented yet");
    }
    public static Cigar leftAlignIndel(Cigar cigar, final byte[] refSeq, final byte[] readSeq, final int refIndex, final int readIndex, final boolean doNotThrowExceptionForMultipleIndels) {
        return leftAlignIndel(cigar, refSeq, readSeq, refIndex, readIndex, 0, doNotThrowExceptionForMultipleIndels);
    }

    private static void ensureLeftAlignmentHasGoodArguments(final Cigar cigar, final byte[] refSeq, final byte[] readSeq, final int refIndex, final int readIndex) {
        Utils.nonNull( cigar );
        Utils.nonNull( refSeq );
        Utils.nonNull( readSeq );
        if ( refIndex < 0 ) throw new IllegalArgumentException("attempting to left align with a reference index less than 0");
        if ( readIndex < 0 ) throw new IllegalArgumentException("attempting to left align with a read index less than 0");
    }

    /**
     * Counts the number of I/D operators
     *
     * @param cigar   cigar to check -- cannot be null
     * @return  non-negative count of indel operators
     */
    private static int countIndelElements(final Cigar cigar) {
        int indelCount = 0;
        for ( CigarElement ce : cigar.getCigarElements() ) {
            if ( ce.getOperator() == CigarOperator.D || ce.getOperator() == CigarOperator.I )
                indelCount++;
        }
        return indelCount;
    }

    /**
     * See the documentation for AlignmentUtils.leftAlignIndel() for more details.
     *
     * This flavor of the left alignment works if and only if the alignment has one - and only one - indel.
     * An exception is thrown if there are no indels or more than 1 indel in the alignment.
     *
     * @param cigar     structure of the original alignment -- cannot be null
     * @param refSeq    reference sequence the read is aligned to
     * @param readSeq   read sequence
     * @param refIndex  0-based alignment start position on ref
     * @param readIndex 0-based alignment start position on read
     * @param leftmostAllowedAlignment left align indel no further left than this index (0-based)
     * @param cleanupCigar if true, we'll cleanup the resulting cigar element, removing 0 length elements and deletions from the first cigar position
     * @return a non-null cigar, in which the single indel is guaranteed to be placed at the leftmost possible position across a repeat (if any)
     */
    public static Cigar leftAlignSingleIndel(Cigar cigar, final byte[] refSeq, final byte[] readSeq, final int refIndex, final int readIndex, final int leftmostAllowedAlignment,final boolean cleanupCigar) {
        ensureLeftAlignmentHasGoodArguments(cigar, refSeq, readSeq, refIndex, readIndex);

        int indexOfIndel = -1;
        for (int i = 0; i < cigar.numCigarElements(); i++) {
            CigarElement ce = cigar.getCigarElement(i);
            if (ce.getOperator() == CigarOperator.D || ce.getOperator() == CigarOperator.I) {
                // if there is more than 1 indel, exception out
                if (indexOfIndel != -1)
                    throw new IllegalArgumentException("attempting to left align a CIGAR that has more than 1 indel in its alignment");
                indexOfIndel = i;
            }
        }

        // if there is no indel, exception out
        if ( indexOfIndel == -1 )
            throw new IllegalArgumentException("attempting to left align a CIGAR that has no indels in its alignment");
        // if the alignment starts with an insertion (so that there is no place on the read to move that insertion further left), we are done
        if ( indexOfIndel == 0 )
            return cigar;

        final int indelLength = cigar.getCigarElement(indexOfIndel).getLength();

        byte[] altString = createIndelString(cigar, indexOfIndel, refSeq, readSeq, refIndex, readIndex);
        if (altString == null)
            return cigar;

        Cigar newCigar = cigar;
        for (int i = 0; i < indelLength; i++) {
            newCigar = moveCigarLeft(newCigar, indexOfIndel);
            if(isIndelAlignedTooFarLeft(newCigar,leftmostAllowedAlignment)) {
                break;
            }
            byte[] newAltString = createIndelString(newCigar, indexOfIndel, refSeq, readSeq, refIndex, readIndex);

            // check to make sure we haven't run off the end of the read
            boolean reachedEndOfRead = cigarHasZeroSizeElement(newCigar);

            if (Arrays.equals(altString, newAltString)) {
                cigar = newCigar;
                i = -1;
                if (reachedEndOfRead)
                    cigar = cleanupCigar ? cleanUpCigar(cigar) : cigar;
            }

            if (reachedEndOfRead)
                break;
        }

        return cigar;
    }
    public static Cigar leftAlignSingleIndel(Cigar cigar, final byte[] refSeq, final byte[] readSeq, final int refIndex, final int readIndex, final boolean cleanupCigar) {
        return leftAlignSingleIndel(cigar, refSeq, readSeq, refIndex, readIndex, 0, cleanupCigar);
    }

    /**
     * Check if cigar aligns indel too far left
     * @param cigar     structure of the original alignment -- cannot be null
     * @param leftmostAllowedAlignment furthest left in cigar indel can be
     * @param
     * @return true is indel is aligned too far left
     */
    protected static boolean isIndelAlignedTooFarLeft(final Cigar cigar,final int leftmostAllowedAlignment) {
        int location=0;
        for (CigarElement element : cigar.getCigarElements() ) {
            if (element.getOperator()==CigarOperator.D || element.getOperator()==CigarOperator.I) {
                return location<leftmostAllowedAlignment;
            }
            if(element.getOperator().consumesReferenceBases()) {
                location += element.getLength();
            }
        }
        return false;
    }


    /**
     * Does one of the elements in cigar have a 0 length?
     *
     * @param c a non-null cigar
     * @return true if any element has 0 size
     */
    protected static boolean cigarHasZeroSizeElement(final Cigar c) {
        for (final CigarElement ce : c.getCigarElements()) {
            if (ce.getLength() == 0)
                return true;
        }
        return false;
    }

    /**
     * Clean up the incoming cigar
     *
     * Removes elements with zero size
     * Clips away beginning deletion operators
     *
     * @param c the cigar string we want to clean up
     * @return a newly allocated, cleaned up Cigar
     */
    public static Cigar cleanUpCigar(final Cigar c) {
        final List<CigarElement> elements = new ArrayList<>(c.numCigarElements() - 1);

        for (final CigarElement ce : c.getCigarElements()) {
            if (ce.getLength() != 0 && (! elements.isEmpty() || ce.getOperator() != CigarOperator.D)) {
                elements.add(ce);
            }
        }

        return new Cigar(elements);
    }

    /**
     * Removing a trailing deletion from the incoming cigar if present
     *
     * @param c the cigar we want to update
     * @return a non-null Cigar
     */
    public static Cigar removeTrailingDeletions(final Cigar c) {

        final List<CigarElement> elements = c.getCigarElements();
        if ( elements.get(elements.size() - 1).getOperator() != CigarOperator.D )
            return c;

        return new Cigar(elements.subList(0, elements.size() - 1));
    }

    /**
     * Move the indel in a given cigar string one base to the left
     *
     * @param cigar          original cigar
     * @param indexOfIndel   the index of the indel cigar element
     * @return non-null cigar with indel moved one base to the left
     */
    private static Cigar moveCigarLeft(Cigar cigar, int indexOfIndel) {
        // get the first few elements
        ArrayList<CigarElement> elements = new ArrayList<>(cigar.numCigarElements());
        for (int i = 0; i < indexOfIndel - 1; i++)
            elements.add(cigar.getCigarElement(i));

        // get the indel element and move it left one base
        CigarElement ce = cigar.getCigarElement(indexOfIndel - 1);
        elements.add(new CigarElement(Math.max(ce.getLength() - 1, 0), ce.getOperator()));
        elements.add(cigar.getCigarElement(indexOfIndel));
        if (indexOfIndel + 1 < cigar.numCigarElements()) {
            ce = cigar.getCigarElement(indexOfIndel + 1);
            elements.add(new CigarElement(ce.getLength() + 1, ce.getOperator()));
        } else {
            elements.add(new CigarElement(1, CigarOperator.M));
        }

        // get the last few elements
        for (int i = indexOfIndel + 2; i < cigar.numCigarElements(); i++)
            elements.add(cigar.getCigarElement(i));
        return new Cigar(elements);
    }

    /**
     * Create the string (really a byte array) representation of an indel-containing cigar against the reference.
     *
     * @param cigar             the indel-containing cigar
     * @param indexOfIndel      the index of the indel cigar element
     * @param refSeq            the reference sequence
     * @param readSeq           the read sequence for the cigar
     * @param refIndex          the starting reference index into refSeq
     * @param readIndex         the starting read index into readSeq
     * @return non-null byte array which is the indel representation against the reference
     */
    private static byte[] createIndelString(final Cigar cigar, final int indexOfIndel, final byte[] refSeq, final byte[] readSeq, int refIndex, int readIndex) {
        CigarElement indel = cigar.getCigarElement(indexOfIndel);
        int indelLength = indel.getLength();

        int totalRefBases = 0;
        for (int i = 0; i < indexOfIndel; i++) {
            CigarElement ce = cigar.getCigarElement(i);
            int length = ce.getLength();

            switch (ce.getOperator()) {
                case M:
                case EQ:
                case X:
                    readIndex += length;
                    refIndex += length;
                    totalRefBases += length;
                    break;
                case S:
                    readIndex += length;
                    break;
                case N:
                    refIndex += length;
                    totalRefBases += length;
                    break;
                default:
                    break;
            }
        }

        // sometimes, when there are very large known indels, we won't have enough reference sequence to cover them
        if (totalRefBases + indelLength > refSeq.length)
            indelLength -= (totalRefBases + indelLength - refSeq.length);

        // the indel-based reference string
        byte[] alt = new byte[refSeq.length + (indelLength * (indel.getOperator() == CigarOperator.D ? -1 : 1))];

        // add the bases before the indel, making sure it's not aligned off the end of the reference
        if (refIndex > alt.length || refIndex > refSeq.length)
            return null;
        System.arraycopy(refSeq, 0, alt, 0, refIndex);
        int currentPos = refIndex;

        // take care of the indel
        if (indel.getOperator() == CigarOperator.D) {
            refIndex += indelLength;
        } else {
            System.arraycopy(readSeq, readIndex, alt, currentPos, indelLength);
            currentPos += indelLength;
        }

        // add the bases after the indel, making sure it's not aligned off the end of the reference
        if (refSeq.length - refIndex > alt.length - currentPos)
            return null;
        System.arraycopy(refSeq, refIndex, alt, currentPos, refSeq.length - refIndex);

        return alt;
    }


    /**
     * Trim cigar down to one that starts at start reference on the left and extends to end on the reference
     *
     * @param cigar a non-null Cigar to trim down
     * @param start Where should we start keeping bases on the reference?  The first position is 0
     * @param end Where should we stop keeping bases on the reference?  The maximum value is cigar.getReferenceLength()
     * @return a new Cigar with reference length == start - end + 1
     */
    public static Cigar trimCigarByReference(final Cigar cigar, final int start, final int end) {
        if ( start < 0 ) throw new IllegalArgumentException("Start must be >= 0 but got " + start);
        if ( end < start ) throw new IllegalArgumentException("End " + end + " is < start start " + start);
        if ( end > cigar.getReferenceLength() ) throw new IllegalArgumentException("End is beyond the cigar's reference length " + end + " for cigar " + cigar );

        final Cigar result = trimCigar(cigar, start, end, true);

        Utils.validate(result.getReferenceLength() == end - start + 1, () -> "trimCigarByReference failure: start " + start + " end " + end + " for " + cigar + " resulted in cigar with wrong size " + result);
        return result;
    }

    /**
     * Trim cigar down to one that starts at start base in the cigar and extends to (inclusive) end base
     *
     * @param cigar a non-null Cigar to trim down
     * @param start Where should we start keeping bases in the cigar?  The first position is 0
     * @param end Where should we stop keeping bases in the cigar?  The maximum value is cigar.getLength()
     * @return a new Cigar containing == start - end + 1 reads
     */
    public static Cigar trimCigarByBases(final Cigar cigar, final int start, final int end) {
        if ( start < 0 ) throw new IllegalArgumentException("Start must be >= 0 but got " + start);
        if ( end < start ) throw new IllegalArgumentException("End " + end + " is < start = " + start);
        if ( end > cigar.getReadLength() ) throw new IllegalArgumentException("End is beyond the cigar's read length " + end + " for cigar " + cigar );

        final Cigar result = trimCigar(cigar, start, end, false);

        final int expectedSize = end - start + 1;
        Utils.validate(result.getReadLength() == expectedSize, () -> "trimCigarByBases failure: start "
                + start + " end " + end + " for " + cigar + " resulted in cigar with wrong size " + result + " with size " + result.getReadLength() + " expected " + expectedSize + " for input cigar " + cigar);
        return result;
    }


    /**
     * Workhorse for trimCigarByBases and trimCigarByReference
     *
     * @param cigar a non-null Cigar to trim down
     * @param start Where should we start keeping bases in the cigar?  The first position is 0
     * @param end Where should we stop keeping bases in the cigar?  The maximum value is cigar.getLength()
     * @param byReference should start and end be intrepreted as position in the reference or the read to trim to/from?
     * @return a non-null cigar
     */
    @SuppressWarnings("fallthrough")
    private static Cigar trimCigar(final Cigar cigar, final int start, final int end, final boolean byReference) {
        final List<CigarElement> newElements = new LinkedList<>();

        int pos = 0;
        for ( final CigarElement elt : cigar.getCigarElements() ) {
            if ( pos > end && (byReference || elt.getOperator() != CigarOperator.D) ) break;

            switch ( elt.getOperator() ) {
                case D:
                    if ( ! byReference ) {
                        if ( pos >= start )
                            newElements.add(elt);
                        break;
                    }
                    // otherwise fall through to the next case
                case EQ: case M: case X:
                    pos = addCigarElements(newElements, pos, start, end, elt);
                    break;
                case S: case I:
                    if ( byReference ) {
                        if ( pos >= start )
                            newElements.add(elt);
                    } else {
                        pos = addCigarElements(newElements, pos, start, end, elt);
                    }
                    break;
                default:
                    throw new IllegalStateException("Cannot handle " + elt);
            }
        }

        return AlignmentUtils.consolidateCigar(new Cigar(newElements));
    }

    /**
     * Helper function for trimCigar that adds cigar elements (of total length X) of elt.op to dest for
     * X bases that fall between start and end, where the last position of the base is pos.
     *
     * The primary use of this function is to create a new cigar element list that contains only
     * elements that occur between start and end bases in an initial cigar.
     *
     * Note that this function may return multiple cigar elements (1M1M etc) that are best consolidated
     * after the fact into a single simpler representation.
     *
     * @param dest we will append our cigar elements to this list
     * @param pos the position (0 indexed) where elt started
     * @param start only include bases that occur >= this position
     * @param end only include bases that occur <= this position
     * @param elt the element we are slicing down
     * @return the position after we've traversed all elt.length bases of elt
     */
    protected static int addCigarElements(final List<CigarElement> dest, int pos, final int start, final int end, final CigarElement elt) {
        final int length = Math.min(pos + elt.getLength() - 1, end) - Math.max(pos, start) + 1;
        if ( length > 0 )
            dest.add(new CigarElement(length, elt.getOperator()));
        return pos + elt.getLength();
    }

    /**
     * Get the offset (base 0) of the first reference aligned base in Cigar that occurs after readStartByBaseOfCigar base of the cigar
     *
     * The main purpose of this routine is to find a good start position for a read given it's cigar.  The real
     * challenge is that the starting base might be inside an insertion, in which case the read actually starts
     * at the next M/EQ/X operator.
     *
     * @param cigar a non-null cigar
     * @param readStartByBaseOfCigar finds the first base after this (0 indexed) that aligns to the reference genome (M, EQ, X)
     * @throws IllegalStateException if no such base can be found
     * @return an offset into cigar
     */
    public static int calcFirstBaseMatchingReferenceInCigar(final Cigar cigar, int readStartByBaseOfCigar) {
        if ( cigar == null ) throw new IllegalArgumentException("cigar cannot be null");
        if ( readStartByBaseOfCigar >= cigar.getReadLength() ) throw new IllegalArgumentException("readStartByBaseOfCigar " + readStartByBaseOfCigar + " must be <= readLength " + cigar.getReadLength());

        int hapOffset = 0, refOffset = 0;
        for ( final CigarElement ce : cigar.getCigarElements() ) {
            for ( int i = 0; i < ce.getLength(); i++ ) {
                switch ( ce.getOperator() ) {
                    case M:case EQ:case X:
                        if ( hapOffset >= readStartByBaseOfCigar )
                            return refOffset;
                        hapOffset++;
                        refOffset++;
                        break;
                    case I: case S:
                        hapOffset++;
                        break;
                    case D:
                        refOffset++;
                        break;
                    default:
                        throw new IllegalStateException("calcFirstBaseMatchingReferenceInCigar does not support cigar " + ce.getOperator() + " in cigar " + cigar);
                }
            }
        }

        throw new IllegalStateException("Never found appropriate matching state for cigar " + cigar + " given start of " + readStartByBaseOfCigar);
    }

    /**
     * Generate a new Cigar that maps the operations of the first cigar through those in a second
     *
     * For example, if first is 5M and the second is 2M1I2M then the result is 2M1I2M.
     * However, if first is 1M2D3M and second is 2M1I3M this results in a cigar X
     *
     * ref   : AC-GTA
     * hap   : ACxGTA  - 2M1I3M
     * read  : A--GTA  - 1M2D3M
     * result: A--GTA => 1M1D3M
     *
     * ref   : ACxG-TA
     * hap   : AC-G-TA  - 2M1D3M
     * read  : AC-GxTA  - 3M1I2M
     * result: AC-GxTA => 2M1D1M1I2M
     *
     * ref   : ACGTA
     * hap   : ACGTA  - 5M
     * read  : A-GTA  - 1M1I3M
     * result: A-GTA => 1M1I3M
     *
     * ref   : ACGTAC
     * hap   : AC---C  - 2M3D1M
     * read  : AC---C  - 3M
     * result: AG---C => 2M3D
     *
     * The constraint here is that both cigars should imply that the result have the same number of
     * reference bases (i.e.g, cigar.getReferenceLength() are equals).
     *
     * @param firstToSecond the cigar mapping hap1 -> hap2
     * @param secondToThird the cigar mapping hap2 -> hap3
     * @return A cigar mapping hap1 -> hap3
     */
    public static Cigar applyCigarToCigar(final Cigar firstToSecond, final Cigar secondToThird) {
        final boolean DEBUG = false;

        final List<CigarElement> newElements = new LinkedList<>();
        final int nElements12 = firstToSecond.numCigarElements();
        final int nElements23 = secondToThird.numCigarElements();

        int cigar12I = 0, cigar23I = 0;
        int elt12I = 0, elt23I = 0;

        while ( cigar12I < nElements12 && cigar23I < nElements23 ) {
            final CigarElement elt12 = firstToSecond.getCigarElement(cigar12I);
            final CigarElement elt23 = secondToThird.getCigarElement(cigar23I);

            final CigarPairTransform transform = getTransformer(elt12.getOperator(), elt23.getOperator());

            if ( DEBUG )
                System.out.printf("Transform %s => %s with elt1 = %d %s @ %d elt2 = %d %s @ %d with transform %s%n",
                        firstToSecond, secondToThird, cigar12I, elt12.getOperator(), elt12I, cigar23I, elt23.getOperator(), elt23I, transform);

            if ( transform.op13 != null ) // skip no ops
                newElements.add(new CigarElement(1, transform.op13));

            elt12I += transform.advance12;
            elt23I += transform.advance23;

            // if have exhausted our current element, advance to the next one
            if ( elt12I == elt12.getLength() ) { cigar12I++; elt12I = 0; }
            if ( elt23I == elt23.getLength() ) { cigar23I++; elt23I = 0; }
        }

        return AlignmentUtils.consolidateCigar(new Cigar(newElements));
    }

    private static CigarPairTransform getTransformer(final CigarOperator op12, final CigarOperator op23) {
        for ( final CigarPairTransform transform : cigarPairTransformers) {
            if ( transform.op12.contains(op12) && transform.op23.contains(op23) )
                return transform;
        }

        throw new IllegalStateException("No transformer for operators " + op12 + " and " + op23);
    }

    /**
     * transformations that project one alignment state through another
     *
     * Think about this as a state machine, where we have:
     *
     * bases3 : xxx A zzz
     * bases2 : xxx B zzz
     * bases1 : xxx C zzz
     *
     * where A, B and C are alignment states of a three way alignment.  We want to capture
     * the transition from operation mapping 1 -> 2 and an operation mapping 2 -> 3 and its
     * associated mapping from 1 -> 3 and the advancement of the cigar states of 1->2 and 2->3.
     *
     * Imagine that A, B, and C are all equivalent (so that op12 = M and op23 = M).  This implies
     * a mapping of 1->3 of M, and in this case the next states to consider in the 3 way alignment
     * are the subsequent states in 1 and 2 (so that advance12 and advance23 are both 1).
     *
     * Obviously not all of the states and their associated transitions are so simple.  Suppose instead
     * that op12 = I, and op23 = M.  What does this look like:
     *
     * bases3 : xxx - A zzz
     * bases2 : xxx - B zzz
     * bases1 : xxx I C zzz
     *
     * It means that op13 must be an insertion (as we have an extra base in 1 thats not present in 2 and
     * so not present in 3).  We advance the cigar in 1 by 1 (as we've consumed one base in 1 for the I)
     * but we haven't yet found the base corresponding to the M of op23.  So we don't advance23.
     */
    private static final class CigarPairTransform {
        private final EnumSet<CigarOperator> op12, op23;
        private final CigarOperator op13;
        private final int advance12, advance23;

        private CigarPairTransform(CigarOperator op12, CigarOperator op23, CigarOperator op13, int advance12, int advance23) {
            this.op12 = getCigarSet(op12);
            this.op23 = getCigarSet(op23);
            this.op13 = op13;
            this.advance12 = advance12;
            this.advance23 = advance23;
        }

        private static EnumSet<CigarOperator> getCigarSet(final CigarOperator masterOp) {
            switch ( masterOp ) {
                case M: return EnumSet.of(CigarOperator.M, CigarOperator.EQ, CigarOperator.X);
                case I: return EnumSet.of(CigarOperator.I, CigarOperator.S);
                case D: return EnumSet.of(CigarOperator.D);
                default: throw new IllegalStateException("Unexpected state " + masterOp);
            }
        }

        @Override
        public String toString() {
            return "CigarPairTransform{" +
                    "op12=" + op12 +
                    ", op23=" + op23 +
                    ", op13=" + op13 +
                    ", advance12=" + advance12 +
                    ", advance23=" + advance23 +
                    '}';
        }
    }


    private static final List<CigarPairTransform> cigarPairTransformers = Arrays.asList(
            //
            // op12 is a match
            //
            // 3: xxx B yyy
            // ^^^^^^^^^^^^
            // 2: xxx M yyy
            // 1: xxx M yyy
            new CigarPairTransform(CigarOperator.M, CigarOperator.M, CigarOperator.M, 1, 1),
            // 3: xxx I yyy
            // ^^^^^^^^^^^^
            // 2: xxx I yyy
            // 1: xxx M yyy
            new CigarPairTransform(CigarOperator.M, CigarOperator.I, CigarOperator.I, 1, 1),
            // 3: xxx D yyy
            // ^^^^^^^^^^^^
            // 2: xxx D yyy
            // 1: xxx M yyy
            new CigarPairTransform(CigarOperator.M, CigarOperator.D, CigarOperator.D, 0, 1),

            //
            // op12 is a deletion
            //
            // 3: xxx D M yyy
            // ^^^^^^^^^^^^
            // 2: xxx M yyy
            // 1: xxx D yyy
            new CigarPairTransform(CigarOperator.D, CigarOperator.M, CigarOperator.D, 1, 1),
            // 3: xxx D2 D1 yyy
            // ^^^^^^^^^^^^
            // 2: xxx D2 yyy
            // 1: xxx D1 yyy
            new CigarPairTransform(CigarOperator.D, CigarOperator.D, CigarOperator.D, 0, 1),
            // 3: xxx X yyy => no-op, we skip emitting anything here
            // ^^^^^^^^^^^^
            // 2: xxx I yyy
            // 1: xxx D yyy
            new CigarPairTransform(CigarOperator.D, CigarOperator.I, null, 1, 1),

            //
            // op12 is a insertion
            //
            // 3: xxx I M yyy
            // ^^^^^^^^^^^^
            // 2: xxx M yyy
            // 1: xxx I yyy
            new CigarPairTransform(CigarOperator.I, CigarOperator.M, CigarOperator.I, 1, 0),
            // 3: xxx I D yyy
            // ^^^^^^^^^^^^
            // 2: xxx D yyy
            // 1: xxx I yyy
            new CigarPairTransform(CigarOperator.I, CigarOperator.D, CigarOperator.I, 1, 0),
            // 3: xxx I1 I2 yyy
            // ^^^^^^^^^^^^
            // 2: xxx I2 yyy
            // 1: xxx I1 yyy
            new CigarPairTransform(CigarOperator.I, CigarOperator.I, CigarOperator.I, 1, 0)
    );
}