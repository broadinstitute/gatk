package org.broadinstitute.hellbender.utils.read;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;

import java.util.*;
import java.util.stream.IntStream;


public final class AlignmentUtils {
    private static final EnumSet<CigarOperator> ALIGNED_TO_GENOME_OPERATORS = EnumSet.of(CigarOperator.M, CigarOperator.EQ, CigarOperator.X);
    private static final EnumSet<CigarOperator> ALIGNED_TO_GENOME_PLUS_SOFTCLIPS = EnumSet.of(CigarOperator.M, CigarOperator.EQ, CigarOperator.X, CigarOperator.S);
    public final static String HAPLOTYPE_TAG = "HC";
    public final static byte GAP_CHARACTER = (byte)'-';

    // cannot be instantiated
    private AlignmentUtils() { }

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

        // compute the smith-waterman alignment of read -> haplotype //TODO use more efficient than the read clipper here
        final GATKRead readMinusSoftClips = ReadClipper.hardClipSoftClippedBases(originalRead);
        final int softClippedBases = originalRead.getLength() - readMinusSoftClips.getLength();
        final SmithWatermanAlignment readToHaplotypeSWAlignment = aligner.align(haplotype.getBases(), readMinusSoftClips.getBases(), CigarUtils.ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS, SWOverhangStrategy.SOFTCLIP);
        if ( readToHaplotypeSWAlignment.getAlignmentOffset() == -1 ) { // START HERE can I get a smith-waterman score?
            // sw can fail (reasons not clear) so if it happens just don't realign the read
            return originalRead;
        }

        final Cigar swCigar = new CigarBuilder().addAll(readToHaplotypeSWAlignment.getCigar()).make();

        // since we're modifying the read we need to clone it
        final GATKRead copiedRead = originalRead.copy();

        // only informative reads are given the haplotype tag to enhance visualization
        if ( isInformative ) {
            copiedRead.setAttribute(HAPLOTYPE_TAG, haplotype.hashCode());
        }

        // compute here the read starts w.r.t. the reference from the SW result and the hap -> ref cigar
        final Cigar rightPaddedHaplotypeVsRefCigar = haplotype.getConsolidatedPaddedCigar(1000);

        // this computes the number of reference bases before the read starts, based on the haplotype vs ref cigar
        // This cigar corresponds exactly to the readToRefCigarRaw, below.  One might wonder whether readToRefCigarRaw and
        // readToRefCigarClean ever imply different starts, which could occur if if the former has a leading deletion.  However,
        // according to the logic of applyCigarToCigar, this can only happen if the read has a leading deletion wrt its best haplotype,
        // which our SW aligner won't do, or if the read starts on a haplotype base that is in a deletion wrt to reference, which is nonsensical
        // since a base that exists is not a deletion.  Thus, there is nothing to worry about, in contrast to below where we do check
        // whether left-alignment shifted the start position.
        final int readStartOnReferenceHaplotype = readStartOnReferenceHaplotype(rightPaddedHaplotypeVsRefCigar, readToHaplotypeSWAlignment.getAlignmentOffset());


        //final int readStartOnReference = referenceStart + haplotype.getAlignmentStartHapwrtRef() + readStartOnHaplotype;

        final int readStartOnReference = referenceStart + haplotype.getAlignmentStartHapwrtRef() + readStartOnReferenceHaplotype;

        // compute the read -> ref alignment by mapping read -> hap -> ref from the
        // SW of read -> hap mapped through the given by hap -> ref

        // this is the sub-cigar of the haplotype-to-ref alignment, with cigar elements before the read start removed.  Elements after the read end are kept.
        final Cigar haplotypeToRef = trimCigarByBases(rightPaddedHaplotypeVsRefCigar, readToHaplotypeSWAlignment.getAlignmentOffset(), rightPaddedHaplotypeVsRefCigar.getReadLength() - 1).getCigar();

        final Cigar readToRefCigar = applyCigarToCigar(swCigar, haplotypeToRef);
        final CigarBuilder.Result leftAlignedReadToRefCigarResult = leftAlignIndels(readToRefCigar, refHaplotype.getBases(), readMinusSoftClips.getBases(), readStartOnReferenceHaplotype);
        final Cigar leftAlignedReadToRefCigar = leftAlignedReadToRefCigarResult.getCigar();
        // it's possible that left-alignment shifted a deletion to the beginning of a read and removed it, shifting the first aligned base to the right
        copiedRead.setPosition(copiedRead.getContig(), readStartOnReference + leftAlignedReadToRefCigarResult.getLeadingDeletionBasesRemoved());

        // the SW Cigar does not contain the hard clips of the original read
        // Here we reconcile the aligned read (that has had any softclips removed) with its softclipped bases
        final Cigar originalCigar = originalRead.getCigar();
        final Cigar newCigar = appendClippedElementsFromCigarToCigar(leftAlignedReadToRefCigar, originalCigar);
        copiedRead.setCigar(newCigar);

        if ( leftAlignedReadToRefCigar.getReadLength() + softClippedBases != copiedRead.getLength() ) {
            throw new GATKException("Cigar " + leftAlignedReadToRefCigar + " with read length " + leftAlignedReadToRefCigar.getReadLength()
                    + " != read length " + copiedRead.getLength() + " for read " + copiedRead.toString() + "\nhapToRef " + haplotypeToRef + " length " + haplotypeToRef.getReadLength() + "/" + haplotypeToRef.getReferenceLength()
                    + "\nreadToHap " + swCigar + " length " + swCigar.getReadLength() + "/" + swCigar.getReferenceLength());
        }
        // assert that the cigar has the same number of elements as the original read

        return copiedRead;
    }

    /**
     * Helper method that handles the work of re-appending clipped bases from the original cigar to the new one
     *
     * @param cigarToHaveClippedElementsAdded cigar to attach softclips to
     * @param originalClippedCigar cigar to check for clipped bases
     * @return a new cigar that has had the clipped elements from the original appended to either end
     */
    @VisibleForTesting
    static Cigar appendClippedElementsFromCigarToCigar(final Cigar cigarToHaveClippedElementsAdded, final Cigar originalClippedCigar) {
        int firstIndex = 0;
        int lastIndex = originalClippedCigar.numCigarElements() - 1;
        CigarElement firstElement = originalClippedCigar.getFirstCigarElement();
        CigarElement lastElement = originalClippedCigar.getLastCigarElement();
        final List<CigarElement> readToRefCigarElementsWithHardClips = new ArrayList<>();
        while (firstElement.getOperator().isClipping() && firstIndex != lastIndex) {
            readToRefCigarElementsWithHardClips.add(firstElement);
            firstElement = originalClippedCigar.getCigarElement(++firstIndex);
        }
        readToRefCigarElementsWithHardClips.addAll(cigarToHaveClippedElementsAdded.getCigarElements());

        final List<CigarElement> endCigarElementsToReverse = new ArrayList<>();
        while (lastElement.getOperator().isClipping() && firstIndex != lastIndex)  {
            endCigarElementsToReverse.add(lastElement);
            lastElement = originalClippedCigar.getCigarElement(--lastIndex);
        }
        // Reverse the order to preserve the original soft/hardclip ordering in mixed clipping cases where softclips precede hardclips
        readToRefCigarElementsWithHardClips.addAll(Lists.reverse(endCigarElementsToReverse));

        return new Cigar(readToRefCigarElementsWithHardClips);
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
     * If the cigar starts with an insertion, the inserted bases are considered as coming before the start position and
     * are therefore excluded from the result.  That is getBasesCoveringRefInterval(0, 3, "ACTTGT", 0, 2I4M) should yield "TTGT".
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

    /**
     * Returns the "IGV View" of all the bases and base qualities in a read aligned to the reference according to the cigar, dropping any bases
     * that might be in the read but aren't in the reference. Any bases that appear in the reference but not the read
     * will be filled in with GAP_CHARACTER values for the read bases and 0's for base qualities to indicate that they don't exist.
     *
     * If the cigar for input read is all matches to the reference then this method will return references to the original
     * read base/base quality byte arrays in the underlying SamRecord in order to save on array allocation/copying performance effects.
     *
     * @param read a read to return aligned to the reference
     * @return A Pair of byte arrays where the left array corresponds to the bases aligned to the reference and right
     *         array corresponds to the baseQualities aligned to the reference.
     */
    public static Pair<byte[], byte[]> getBasesAndBaseQualitiesAlignedOneToOne(final GATKRead read) {
        return getBasesAndBaseQualitiesAlignedOneToOne(read, GAP_CHARACTER, (byte)0);
    }

    private static Pair<byte[], byte[]> getBasesAndBaseQualitiesAlignedOneToOne(final GATKRead read, final byte basePadCharacter, final byte qualityPadCharacter) {
        Utils.nonNull(read);
        // As this code is performance sensitive in the HaplotypeCaller, we elect to use the noCopy versions of these getters.
        // We can do this because we don't mutate base or quality arrays in this method or in its accessors
        final byte[] bases = read.getBasesNoCopy();
        final byte[] baseQualities = read.getBaseQualitiesNoCopy();
        final int numCigarElements = read.numCigarElements();
        boolean sawIndel = false;

        // Check if the cigar contains indels
        // Note that we don't call ContainsOperator() here twice to avoid the performance hit of building stream iterators twice
        for (int i = 0; i < numCigarElements; i++) {
            final CigarOperator e = read.getCigarElement(i).getOperator();
            if (e == CigarOperator.INSERTION || e == CigarOperator.DELETION) {
                sawIndel = true;
                break;
            }
        }
        if (!sawIndel) {
            return new ImmutablePair<>(bases, baseQualities);
        }
        else {
            final int numberRefBasesIncludingSoftclips = CigarUtils.countRefBasesAndSoftClips(read.getCigarElements(), 0, numCigarElements);
            final byte[] paddedBases = new byte[numberRefBasesIncludingSoftclips];
            final byte[] paddedBaseQualities = new byte[numberRefBasesIncludingSoftclips];
            int literalPos = 0;
            int paddedPos = 0;
            for ( int i = 0; i < numCigarElements; i++ ) {
                final CigarElement ce = read.getCigarElement(i);
                final CigarOperator co = ce.getOperator();
                if (co.consumesReadBases()) {
                    if (!co.consumesReferenceBases()) {
                        literalPos += ce.getLength();  //skip inserted bases
                    }
                    else {
                        System.arraycopy(bases, literalPos, paddedBases, paddedPos, ce.getLength());
                        System.arraycopy(baseQualities, literalPos, paddedBaseQualities, paddedPos, ce.getLength());
                        literalPos += ce.getLength();
                        paddedPos += ce.getLength();
                    }
                }
                else if (co.consumesReferenceBases()) {
                    for ( int j = 0; j < ce.getLength(); j++ ) {  //pad deleted bases
                        paddedBases[paddedPos] = basePadCharacter;
                        paddedBaseQualities[paddedPos] = qualityPadCharacter;
                        paddedPos++;
                    }
                }
            }
            return new ImmutablePair<>(paddedBases, paddedBaseQualities);
        }
    }

    /**
     * Returns the number of bases in a read minus the number of softclipped bases.
     */
    public static int unclippedReadLength(final GATKRead read) {
        int softClippedBases = 0;
        for (CigarElement element : read.getCigarElements()) {
            if (element.getOperator()== CigarOperator.SOFT_CLIP) {
                softClippedBases+= element.getLength();
            }
        }
        return read.getLength() - softClippedBases;
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
    public static int countHighQualitySoftClips(final GATKRead read, final byte qualThreshold ) {
        Utils.nonNull(read);
        ParamUtils.isPositiveOrZero(qualThreshold, "Expected qualThreshold to be positive");

        final Cigar cigar = read.getCigar();
        if ( cigar == null ) {  // the read is unmapped
            return 0;
        }

        int numHQSoftClips = 0;
        int alignPos = 0;
        for ( final CigarElement elem : read.getCigarElements() ) {
            final int elementLength = elem.getLength();
            final CigarOperator operator = elem.getOperator();

            if (operator == CigarOperator.SOFT_CLIP) {
                for (int n = 0; n < elementLength; n++) {
                    if( read.getBaseQuality(alignPos++) > qualThreshold ) {
                        numHQSoftClips++;
                    }
                }
            } else if (operator.consumesReadBases()) {
                alignPos += elementLength;
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

    private static int lengthOnRead(final CigarElement element) {
        return element.getOperator().consumesReadBases() ? element.getLength() : 0;
    }

    private static int lengthOnReference(final CigarElement element) {
        return element.getOperator().consumesReferenceBases() ? element.getLength() : 0;
    }

    /**
     * Takes the alignment of the read sequence <code>readSeq</code> to the reference sequence <code>refSeq</code>
     * starting at 0-based position <code>readStart</code> on the <code>ref</code> and specified by its <code>cigar</code>.
     * <p/>
     * If the alignment has one or more indels, this method attempts to move them left across a stretch of repetitive bases.
     * For instance, if the original cigar specifies that (any) one AT is deleted from a repeat sequence TATATATA, the output
     * cigar will always mark the leftmost AT as deleted. If there is no indel in the original cigar or if the indel position
     * is determined unambiguously (i.e. inserted/deleted sequence is not repeated), the original cigar is returned.
     *
     * Soft-clipped bases in the cigar are presumed to correspond to bases in the byte[] of read sequence.  That is, this method
     * assumes that inputs are precise about the distinction between hard clips (removed from the read sequence) and soft clips
     * (kept in the read sequence but not aligned).  For example, with the inputs {cigar: 2S2M2I, read sequence: TTAAAA, ref sequence: GGAA, read start: 2}
     * the method lines up the AAAA (2M2I) of the read with the AA of the ref and left-aligns the indel to yield a cigar of
     * 2S2I2M.
     *
     * @param cigar     structure of the original alignment
     * @param ref    reference sequence the read is aligned to
     * @param read   read sequence
     * @param readStart  0-based position on ref of the first aligned base in the read sequence
     * @return a non-null cigar, in which the indels are guaranteed to be placed at the leftmost possible position across a repeat (if any)
     */
    public static CigarBuilder.Result leftAlignIndels(final Cigar cigar, final byte[] ref, final byte[] read, final int readStart) {
        ParamUtils.isPositiveOrZero(readStart, "read start within reference base array must be non-negative");

        if (cigar.getCigarElements().stream().noneMatch(elem -> elem.getOperator().isIndel())) {
            return new CigarBuilder.Result(cigar, 0, 0);
        }

        // we need reference bases from the start of the read to the rightmost indel
        final int lastIndel = IntStream.range(0, cigar.numCigarElements()).filter(n -> cigar.getCigarElement(n).getOperator().isIndel()).max().getAsInt();
        final int necessaryRefLength = readStart + cigar.getCigarElements().stream().limit(lastIndel + 1).mapToInt(e -> lengthOnReference(e)).sum();
        Utils.validateArg(necessaryRefLength <= ref.length, "read goes past end of reference");

        // at this point, we are one base past the end of the read.  Now we traverse the cigar from right to left
        final List<CigarElement> resultRightToLeft = new ArrayList<>();
        final int refLength = cigar.getReferenceLength();
        final IndexRange refIndelRange = new IndexRange(readStart + refLength, readStart + refLength);
        final IndexRange readIndelRange = new IndexRange(read.length,read.length);
        for (int n = cigar.numCigarElements() - 1; n >= 0; n--) {
            final CigarElement element = cigar.getCigarElement(n);
            // if it's an indel, just accumulate the read and ref bases consumed.  We won't shift the indel until we hit an alignment
            // block or the read start.
            if (element.getOperator().isIndel()) {
                refIndelRange.shiftStartLeft(lengthOnReference(element));
                readIndelRange.shiftStartLeft(lengthOnRead(element));
            } else if (refIndelRange.size() == 0 && readIndelRange.size() == 0) {   // no indel, just add the cigar element to the result
                resultRightToLeft.add(element);
                refIndelRange.shiftLeft(lengthOnReference(element));
                readIndelRange.shiftLeft(lengthOnRead(element));
            } else {    // there's an indel that we need to left-align
                // we can left-align into match cigar elements but not into clips
                final int maxShift = element.getOperator().isAlignment() ? element.getLength() : 0;
                final Pair<Integer, Integer> shifts = normalizeAlleles(Arrays.asList(ref, read), Arrays.asList(refIndelRange, readIndelRange), maxShift, true);

                // account for new match alignments on the right due to left-alignment
                resultRightToLeft.add(new CigarElement(shifts.getRight(), CigarOperator.MATCH_OR_MISMATCH));

                // emit if we didn't go all the way to the start of an alignment block OR we have reached clips OR we have reached the start of the cigar
                final boolean emitIndel = n == 0 || shifts.getLeft() < maxShift || !element.getOperator().isAlignment();
                final int newMatchOnLeftDueToTrimming = shifts.getLeft() < 0 ? -shifts.getLeft() : 0;   // we may have actually shifted right to make the alleles parsimonious
                final int remainingBasesOnLeft = shifts.getLeft() < 0 ? element.getLength() : (element.getLength() - shifts.getLeft());

                if (emitIndel) {  // some of this alignment block remains after left-alignment -- emit the indel
                    resultRightToLeft.add(new CigarElement(refIndelRange.size(), CigarOperator.DELETION));
                    resultRightToLeft.add(new CigarElement(readIndelRange.size(), CigarOperator.INSERTION));
                    refIndelRange.shiftEndLeft(refIndelRange.size());       // ref indel range is now empty and points to start of left-aligned indel
                    readIndelRange.shiftEndLeft(readIndelRange.size());     // read indel range is now empty and points to start of left-aligned indel

                    refIndelRange.shiftLeft(newMatchOnLeftDueToTrimming + (element.getOperator().consumesReferenceBases() ? remainingBasesOnLeft : 0));
                    readIndelRange.shiftLeft(newMatchOnLeftDueToTrimming + (element.getOperator().consumesReadBases() ? remainingBasesOnLeft : 0));
                    // now read and ref indel ranges are empty and point to end of element preceding this block
                }
                resultRightToLeft.add(new CigarElement(newMatchOnLeftDueToTrimming, CigarOperator.MATCH_OR_MISMATCH));
                resultRightToLeft.add(new CigarElement(remainingBasesOnLeft, element.getOperator()));
            }
        }

        // account for any indels at the start of the cigar that weren't processed because they have no adjacent non-indel element to the left
        resultRightToLeft.add(new CigarElement(refIndelRange.size(), CigarOperator.DELETION));
        resultRightToLeft.add(new CigarElement(readIndelRange.size(), CigarOperator.INSERTION));

        Utils.validateArg(readIndelRange.getStart() == 0, "Given cigar does not account for all bases of the read");
        return new CigarBuilder().addAll(Lists.reverse(resultRightToLeft)).makeAndRecordDeletionsRemovedResult();
    }

    /**
     *  Example usage:  reference = GAAT, read = GAAAT (insertion of one A) and we initially consider the insertion of the A to occur before
     *  the T.  Thus the reference range of this allele is [3,3) (no bases) and the read range is [3,4).  This will be left-aligned so that
     *  the insertion occurs after the G, so that the ranges become [1,1) and [1,2) and the returned shifts are 2 bases for both the start and end
     *  of the range.
     *
     *  If the given allele ranges are not parsimonious, for example [3,4) and [3,5) in the above example to include the common T in both alleles,
     *  the resulting ranges will be shifted by different amounts.  In this case, the shifts are 2 bases in the front and 3 at the end.
     *
     *  Note that we use the convention that the ref allele in the case of an alt insertion, or the alt allele in case of a deletion, is represented
     *  by [n,n) where n is the last aligned coordinate before the indel.  This makes sense when you think in terms of alignment CIGARs:
     *  suppose for example we have a 5M5I5M read with start 100.  Then the match bases are [100,105) on the ref and [0,5) on the read and the inserted bases are
     *  [5,10) on the read and [5,5) on the reference.
     *
     * @param sequences bases of sequences containing different alleles -- could be reference, a haplotype, a read, or subsequences thereof
     * @param bounds    initial ranges (inclusive start, exclusive end) of alleles in same order as {@code sequences}
     *                  Note that this method adjusts these ranges as a side effect
     * @param maxShift  maximum allowable shift left.  This may be less than the amount demanded by the array bounds.  For example, when
     *                  left-aligning a read with multiple indels, we don't want to realign one indel past another (if they "collide" we merge
     *                  them into a single indel and continue -- see {@link AlignmentUtils::leftAlignIndels}
     * @param trim      If true, remove redundant shared bases at the start and end of alleles
     * @return          The number of bases the alleles were shifted left such that they still represented the same event.
     */
    public static Pair<Integer, Integer> normalizeAlleles(final List<byte[]> sequences, final List<IndexRange> bounds, final int maxShift, final boolean trim) {
        Utils.nonEmpty(sequences);
        Utils.validateArg(sequences.size() == bounds.size(), "Must have one initial allele range per sequence");
        bounds.forEach(bound -> Utils.validateArg(maxShift <= bound.getStart(), "maxShift goes past the start of a sequence"));

        int startShift = 0;
        int endShift = 0;

        // consume any redundant shared bases at the end of the alleles
        int minSize = bounds.stream().mapToInt(IndexRange::size).min().getAsInt();
        while (trim && minSize > 0 && lastBaseOnRightIsSame(sequences, bounds)) {
            bounds.forEach(bound -> bound.shiftEndLeft(1));
            minSize--;
            endShift++;
        }

        while (trim && minSize > 0 && firstBaseOnLeftIsSame(sequences, bounds)) {
            bounds.forEach(bound -> bound.shiftStart(1));
            minSize--;
            startShift--;
        }

        // we shift left as long as the last bases on the right are equal among all sequences and the next bases on the left are all equal.
        // if a sequence is empty (eg the reference relative to an insertion alt allele) the last base on the right is the next base on the left
        while( startShift < maxShift && nextBaseOnLeftIsSame(sequences, bounds) && lastBaseOnRightIsSame(sequences, bounds)) {
                bounds.forEach(b -> b.shiftLeft(1));
                startShift++;
                endShift++;
        }

        return ImmutablePair.of(startShift, endShift);
    }

    // do all sequences share a common base at the end of the given index range
    private static boolean lastBaseOnRightIsSame(List<byte[]> sequences, List<IndexRange> bounds) {
        final byte lastBaseOnRight = sequences.get(0)[bounds.get(0).getEnd() - 1];
        for(int n = 0; n < sequences.size(); n++) {
            if (sequences.get(n)[bounds.get(n).getEnd() - 1] != lastBaseOnRight) {
                return false;
            }
        }
        return true;
    }

    // do all sequences share a common first base
    private static boolean firstBaseOnLeftIsSame(final List<byte[]> sequences, final List<IndexRange> bounds) {
        final byte nextBaseOnLeft = sequences.get(0)[bounds.get(0).getStart()];
        for(int n = 0; n < sequences.size(); n++) {
            if (sequences.get(n)[bounds.get(n).getStart()] != nextBaseOnLeft) {
                return false;
            }
        }
        return true;
    }

    // do all sequences share a common base just before the given index ranges
    private static boolean nextBaseOnLeftIsSame(final List<byte[]> sequences, final List<IndexRange> bounds) {
        final byte nextBaseOnLeft = sequences.get(0)[bounds.get(0).getStart() - 1];
        for(int n = 0; n < sequences.size(); n++) {
            if (sequences.get(n)[bounds.get(n).getStart() - 1] != nextBaseOnLeft) {
                return false;
            }
        }
        return true;
    }

    /**
     * Given a read's first aligned base on an alt haplotype, find the first aligned base in the reference haplotype.  This
     * method assumes that the alt haplotype and reference haplotype start at the same place.  That is, the alt haplotype
     * starts at index 0 within the reference base array.
     *
     * @param haplotypeVsRefCigar
     * @param readStartOnHaplotype
     * @return
     */
    @VisibleForTesting
    static int readStartOnReferenceHaplotype(final Cigar haplotypeVsRefCigar, final int readStartOnHaplotype) {
        if (readStartOnHaplotype == 0) {
            return 0;
        }
        // move forward in the haplotype vs ref cigar until we have consumed enough haplotype bases to reach the read start
        // the number of reference bases consumed during this traversal gives us the reference start
        int refBasesConsumedBeforeReadStart = 0;
        int haplotypeBasesConsumed = 0;
        for (final CigarElement element : haplotypeVsRefCigar) {
            refBasesConsumedBeforeReadStart += lengthOnReference(element);
            haplotypeBasesConsumed += lengthOnRead(element);

            if (haplotypeBasesConsumed >= readStartOnHaplotype) {
                final int excess = element.getOperator().consumesReferenceBases() ? haplotypeBasesConsumed - readStartOnHaplotype : 0;
                return refBasesConsumedBeforeReadStart - excess;
            }
        }

        throw new IllegalArgumentException("Cigar doesn't reach the read start");
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
     * Does cigar start or end with a deletion operation?
     *
     * WARNING: there is usually no good reason to use this method, because one should handle the leading and trailing
     * deletion via the {@link CigarBuilder} class.  Do not use this method when you can instead use {@link CigarBuilder}.
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
     * Trim cigar down to one that starts at start reference on the left and extends to end on the reference
     *
     * @param cigar a non-null Cigar to trim down
     * @param start Where should we start keeping bases on the reference?  The first position is 0
     * @param end Where should we stop keeping bases on the reference?  The maximum value is cigar.getReferenceLength()
     * @return a new Cigar with reference length == start - end + 1
     */
    public static CigarBuilder.Result trimCigarByReference(final Cigar cigar, final int start, final int end) {
        return trimCigar(cigar, start, end, true);
    }

    /**
     * Trim cigar down to one that starts at start base in the cigar and extends to (inclusive) end base
     *
     * @param cigar a non-null Cigar to trim down
     * @param start Where should we start keeping bases in the cigar (inclusive)?  The first position is 0
     * @param end Where should we stop keeping bases in the cigar (inclusive)?  The maximum value is cigar.getLength() - 1
     * @return a new Cigar containing == start - end + 1 reads
     */
    public static CigarBuilder.Result trimCigarByBases(final Cigar cigar, final int start, final int end) {
        return trimCigar(cigar, start, end, false);
    }


    /**
     * Workhorse for trimCigarByBases and trimCigarByReference
     *
     * @param cigar a non-null Cigar to trim down
     * @param start Where should we start keeping bases in the cigar (inclusive)?  The first position is 0
     * @param end Where should we stop keeping bases in the cigar (inclusive)?  The maximum value is cigar.getLength() - 1
     * @param byReference should start and end be interpreted as position in the reference or the read to trim to/from?
     * @return a non-null cigar
     */
    @SuppressWarnings("fallthrough")
    private static CigarBuilder.Result trimCigar(final Cigar cigar, final int start, final int end, final boolean byReference) {
        ParamUtils.isPositiveOrZero(start, "start position can't be negative");
        Utils.validateArg(end >= start, () -> "end " + end + " is before start " + start);
        final CigarBuilder newElements = new CigarBuilder();

        // these variables track the inclusive start and exclusive end of the current cigar element in reference (if byReference) or read (otherwise) coordinates
        int elementStart;           // inclusive
        int elementEnd = 0;         // exclusive -- start of next element
        for ( final CigarElement elt : cigar.getCigarElements() ) {
            elementStart = elementEnd;
            elementEnd = elementStart + (byReference ? lengthOnReference(elt) : lengthOnRead(elt));

            // we are careful to include zero-length elements at both ends, that is, elements with elementStart == elementEnd == start and elementStart == elementEnd == end + 1
            if (elementEnd < start || (elementEnd == start && elementStart < start)) {
                continue;
            } else if (elementStart > end && elementEnd > end + 1) {
                break;
            }

            final int overlapLength = elementEnd == elementStart ? elt.getLength() : Math.min(end + 1, elementEnd) - Math.max(start, elementStart);
            newElements.add(new CigarElement(overlapLength, elt.getOperator()));
        }
        Utils.validateArg(elementEnd > end, () -> "cigar elements don't reach end position (inclusive) " + end);

        return newElements.makeAndRecordDeletionsRemovedResult();
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
        final CigarBuilder newElements = new CigarBuilder();
        final int nElements12 = firstToSecond.numCigarElements();
        final int nElements23 = secondToThird.numCigarElements();

        int cigar12I = 0, cigar23I = 0;
        int elt12I = 0, elt23I = 0;

        while ( cigar12I < nElements12 && cigar23I < nElements23 ) {
            final CigarElement elt12 = firstToSecond.getCigarElement(cigar12I);
            final CigarElement elt23 = secondToThird.getCigarElement(cigar23I);

            final CigarPairTransform transform = getTransformer(elt12.getOperator(), elt23.getOperator());

            if ( transform.op13 != null ) // skip no ops
                newElements.add(new CigarElement(1, transform.op13));

            elt12I += transform.advance12;
            elt23I += transform.advance23;

            // if have exhausted our current element, advance to the next one
            if ( elt12I == elt12.getLength() ) { cigar12I++; elt12I = 0; }
            if ( elt23I == elt23.getLength() ) { cigar23I++; elt23I = 0; }
        }

        return newElements.make();
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