package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.*;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;

import java.util.*;

public final class ArtificialReadUtils {
    public static final int DEFAULT_READ_LENGTH = 50;

    private static final String DEFAULT_READ_GROUP_PREFIX = "ReadGroup";
    private static final String DEFAULT_PLATFORM_UNIT_PREFIX = "Lane";
    private static final String DEFAULT_PLATFORM_PREFIX = "Platform";
    private static final String DEFAULT_SAMPLE_NAME = "SampleX";

    /**
     * Creates an artificial sam header, matching the parameters, chromosomes which will be labeled chr1, chr2, etc
     *
     * @param numberOfChromosomes the number of chromosomes to create
     * @param startingChromosome  the starting number for the chromosome (most likely set to 1)
     * @param chromosomeSize      the length of each chromosome
     * @return
     */
    public static SAMFileHeader createArtificialSamHeader(int numberOfChromosomes, int startingChromosome, int chromosomeSize) {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        SAMSequenceDictionary dict = new SAMSequenceDictionary();
        // make up some sequence records
        for (int x = startingChromosome; x < startingChromosome + numberOfChromosomes; x++) {
            SAMSequenceRecord rec = new SAMSequenceRecord( Integer.toString(x), chromosomeSize /* size */);
            rec.setSequenceLength(chromosomeSize);
            dict.addSequence(rec);
        }
        header.setSequenceDictionary(dict);
        return header;
    }

    /**
     * Creates an artificial sam header, matching the parameters, chromosomes which will be labeled chr1, chr2, etc
     * It also adds read groups.
     *
     * @param numberOfChromosomes the number of chromosomes to create
     * @param startingChromosome  the starting number for the chromosome (most likely set to 1)
     * @param chromosomeSize      the length of each chromosome
     * @param groupCount          the number of groups to make
     */
    public static SAMFileHeader createArtificialSamHeaderWithGroups(int numberOfChromosomes, int startingChromosome, int chromosomeSize, int groupCount) {
        final SAMFileHeader header = createArtificialSamHeader(numberOfChromosomes, startingChromosome, chromosomeSize);

        final List<SAMReadGroupRecord> readGroups = new ArrayList<>();
        for (int i = 0; i < groupCount; i++) {
            SAMReadGroupRecord rec = new SAMReadGroupRecord(DEFAULT_READ_GROUP_PREFIX + i);
            rec.setSample(DEFAULT_SAMPLE_NAME);
            readGroups.add(rec);
        }
        header.setReadGroups(readGroups);

        for (int i = 0; i < groupCount; i++) {
            final SAMReadGroupRecord groupRecord = header.getReadGroup(readGroups.get(i).getId());
            groupRecord.setPlatform(DEFAULT_PLATFORM_PREFIX + ((i % 2)+1));
            groupRecord.setPlatformUnit(DEFAULT_PLATFORM_UNIT_PREFIX + ((i % 3)+1));
        }
        return header;
    }

    /**
     * Creates an artificial sam header with standard test parameters
     *
     * @return the sam header
     */
    public static SAMFileHeader createArtificialSamHeader() {
        return createArtificialSamHeader(22, 1, 1000000);
    }

    public static SAMFileHeader createArtificialSamHeaderWithReadGroup( final SAMReadGroupRecord readGroup ) {
        final SAMFileHeader header = createArtificialSamHeader();
        header.addReadGroup(readGroup);
        return header;
    }

    /**
     * Create an artificial GATKRead based on the parameters.  The cigar string will be *M, where * is the length of the read
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param refIndex       the reference index, i.e. what chromosome to associate it with
     * @param alignmentStart where to start the alignment
     * @param length         the length of the read
     * @return the artificial GATKRead
     */
    public static GATKRead createArtificialRead(SAMFileHeader header, String name, int refIndex, int alignmentStart, int length) {
        return new SAMRecordToGATKReadAdapter(createArtificialSAMRecord(header, name, refIndex, alignmentStart, length));
    }

    /**
     * Create an artificial GATKRead based on the parameters.  The cigar string will be *M, where * is the length of the read
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param refIndex       the reference index, i.e. what chromosome to associate it with
     * @param alignmentStart where to start the alignment
     * @param bases          the sequence of the read
     * @param qual           the qualities of the read
     * @return the artificial GATKRead
     */
    public static GATKRead createArtificialRead(SAMFileHeader header, String name, int refIndex, int alignmentStart, byte[] bases, byte[] qual) {
        return new SAMRecordToGATKReadAdapter(createArtificialSAMRecord(header, name, refIndex, alignmentStart, bases, qual));
    }

    /**
     * Create an artificial GATKRead based on the parameters.  The cigar string will be *M, where * is the length of the read
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param contig         the contig to which the read is aligned
     * @param alignmentStart where to start the alignment
     * @param bases          the sequence of the read
     * @param qual           the qualities of the read
     * @return the artificial GATKRead
     */
    public static GATKRead createArtificialRead(SAMFileHeader header, String name, String contig, int alignmentStart, byte[] bases, byte[] qual) {
        return new SAMRecordToGATKReadAdapter(createArtificialSAMRecord(header, name, header.getSequenceIndex(contig), alignmentStart, bases, qual));
    }

    /**
     * Create an artificial GATKRead based on the parameters
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param refIndex       the reference index, i.e. what chromosome to associate it with
     * @param alignmentStart where to start the alignment
     * @param bases          the sequence of the read
     * @param qual           the qualities of the read
     * @param cigar          the cigar string of the read
     * @return the artificial GATKRead
     */
    public static GATKRead createArtificialRead(SAMFileHeader header, String name, int refIndex, int alignmentStart, byte[] bases, byte[] qual, String cigar) {
        return new SAMRecordToGATKReadAdapter(createArtificialSAMRecord(header, name, refIndex, alignmentStart, bases, qual, cigar));
    }


    /**
     * Create an artificial GATKRead based on the parameters
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param contig         the contig to which the read is aligned
     * @param alignmentStart where to start the alignment
     * @param bases          the sequence of the read
     * @param qual           the qualities of the read
     * @param cigar          the cigar string of the read
     * @return the artificial GATKRead
     */
    public static GATKRead createArtificialRead(SAMFileHeader header, String name, String contig, int alignmentStart, byte[] bases, byte[] qual, String cigar) {
        return new SAMRecordToGATKReadAdapter(createArtificialSAMRecord(header, name, header.getSequenceIndex(contig), alignmentStart, bases, qual, cigar));
    }

    /**
     * Create an artificial GATKRead with the following default parameters :
     * header:
     * numberOfChromosomes = 1
     * startingChromosome = 1
     * chromosomeSize = 1000000
     * read:
     * name = "default_read"
     * refIndex = 0
     * alignmentStart = 10000
     *
     * @param header SAM header for the read
     * @param bases the sequence of the read
     * @param qual  the qualities of the read
     * @param cigar the cigar string of the read
     * @return the artificial GATKRead
     */
    public static GATKRead createArtificialRead(final SAMFileHeader header, final byte[] bases, final byte[] qual, final String cigar) {
        return new SAMRecordToGATKReadAdapter(createArtificialSAMRecord(header, "default_read", 0, 10000, bases, qual, cigar));
    }

    public static GATKRead createArtificialRead(final byte[] bases, final byte[] qual, final String cigar) {
        SAMFileHeader header = createArtificialSamHeader();
        return new SAMRecordToGATKReadAdapter(createArtificialSAMRecord(header, "default_read", 0, 10000, bases, qual, cigar));
    }

    public static GATKRead createArtificialRead(final SAMFileHeader header, final String cigarString) {
        return new SAMRecordToGATKReadAdapter(createArtificialSAMRecord(header, TextCigarCodec.decode(cigarString)));
    }

    public static GATKRead createArtificialRead(final SAMFileHeader header, final Cigar cigar) {
        return new SAMRecordToGATKReadAdapter(createArtificialSAMRecord(header, cigar));
    }

    public static GATKRead createArtificialRead(final Cigar cigar) {
        return new SAMRecordToGATKReadAdapter(createArtificialSAMRecord(cigar));
    }

    /**
     * Makes a new read with a name that is unique (so that it will return false to equals(otherRead)
     */
    public static GATKRead createUniqueArtificialRead(final Cigar cigar) {
        return new SAMRecordToGATKReadAdapter(createUniqueArtificialSAMRecord(cigar));
    }

    public static GATKRead createArtificialRead(final Cigar cigar, final String name) {
        return new SAMRecordToGATKReadAdapter(createArtificialSAMRecord(createArtificialSamHeader(), cigar, name));
    }

    public static GATKRead createArtificialRead(final String cigarString) {
        return new SAMRecordToGATKReadAdapter(createArtificialSAMRecord(TextCigarCodec.decode(cigarString)));
    }

    public static GATKRead createArtificialUnmappedRead(final SAMFileHeader header, final byte[] bases, final byte[] qual) {
        final SAMRecord read = new SAMRecord(header);
        read.setReadUnmappedFlag(true);
        read.setReadBases(bases);
        read.setBaseQualities(qual);

        return new SAMRecordToGATKReadAdapter(read);
    }

    public static GATKRead createArtificialUnmappedReadWithAssignedPosition(final SAMFileHeader header, final String contig, final int alignmentStart, final byte[] bases, final byte[] qual) {
        final SAMRecord read = new SAMRecord(header);
        read.setReferenceName(contig);
        read.setAlignmentStart(alignmentStart);
        read.setReadUnmappedFlag(true);
        read.setReadBases(bases);
        read.setBaseQualities(qual);

        return new SAMRecordToGATKReadAdapter(read);
    }

    /**
     * Makes a new read with a name that is unique (so that it will return false to equals(otherRead)
     */
    public static GATKRead createUniqueArtificialRead(final String cigarString) {
        return new SAMRecordToGATKReadAdapter(createUniqueArtificialSAMRecord(TextCigarCodec.decode(cigarString)));
    }

    /**
     * Creates an artificial GATKRead backed by a SAMRecord.
     *
     * The read will consist of the specified number of Q30 'A' bases, and will be
     * mapped to contig "1" at the specified start position.
     *
     * @param name name of the new read
     * @param start start position of the new read
     * @param length number of bases in the new read
     * @return an artificial GATKRead backed by a SAMRecord.
     */
    public static GATKRead createSamBackedRead( final String name, final int start, final int length ) {
        return createSamBackedRead(name, "1", start, length);
    }

    /**
     * Creates an artificial GATKRead backed by a SAMRecord.
     *
     * The read will consist of the specified number of Q30 'A' bases, and will be
     * mapped to the specified contig at the specified start position.
     *
     * @param name name of the new read
     * @param contig contig the new read is mapped to
     * @param start start position of the new read
     * @param length number of bases in the new read
     * @return an artificial GATKRead backed by a SAMRecord.
     */
    public static GATKRead createSamBackedRead( final String name, final String contig, final int start, final int length ) {
        final SAMFileHeader header = createArtificialSamHeader();
        final byte[] bases = Utils.dupBytes((byte)'A', length);
        final byte[] quals = Utils.dupBytes((byte) 30, length);

        final SAMRecord sam = createArtificialSAMRecord(header, bases, quals, length + "M");
        sam.setReadName(name);
        sam.setReferenceName(contig);
        sam.setAlignmentStart(start);
        return new SAMRecordToGATKReadAdapter(sam);
    }

    /**
     * Creates an artificial GATKRead backed by a SAMRecord with no header.
     *
     * The read will consist of the specified number of Q30 'A' bases, and will be
     * mapped to the specified contig at the specified start position.
     *
     * @param name name of the new read
     * @param contig contig the new read is mapped to
     * @param start start position of the new read
     * @param length number of bases in the new read
     * @return an artificial GATKRead backed by a SAMRecord.
     */
    public static GATKRead createHeaderlessSamBackedRead( final String name, final String contig, final int start, final int length ) {
        GATKRead read = createSamBackedRead(name, contig, start, length);
        ((SAMRecordToGATKReadAdapter) read).getEncapsulatedSamRecord().setHeaderStrict(null);
        return read;
    }

    /**
     * Create an artificial SAMRecord based on the parameters.  The cigar string will be *M, where * is the length of the read
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param refIndex       the reference index, i.e. what chromosome to associate it with
     * @param alignmentStart where to start the alignment
     * @param length         the length of the read
     * @return the artificial SAMRecord
     */
    public static SAMRecord createArtificialSAMRecord(SAMFileHeader header, String name, int refIndex, int alignmentStart, int length) {
        if ((refIndex == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && alignmentStart != SAMRecord.NO_ALIGNMENT_START) ||
                (refIndex != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && alignmentStart == SAMRecord.NO_ALIGNMENT_START))
            throw new IllegalArgumentException("Invalid alignment start for artificial read, start = " + alignmentStart);
        SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceIndex(refIndex);
        record.setAlignmentStart(alignmentStart);
        List<CigarElement> elements = new ArrayList<>();
        elements.add(new CigarElement(length, CigarOperator.characterToEnum('M')));
        record.setCigar(new Cigar(elements));
        record.setProperPairFlag(false);

        // our reads and quals are all 'A's by default
        byte[] c = new byte[length];
        byte[] q = new byte[length];
        for (int x = 0; x < length; x++)
            c[x] = q[x] = 'A';
        record.setReadBases(c);
        record.setBaseQualities(q);

        if (refIndex == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
            record.setReadUnmappedFlag(true);
        }

        return record;
    }

    /**
     * Create an artificial SAMRecord based on the parameters.  The cigar string will be *M, where * is the length of the read
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param refIndex       the reference index, i.e. what chromosome to associate it with
     * @param alignmentStart where to start the alignment
     * @param bases          the sequence of the read
     * @param qual           the qualities of the read
     * @return the artificial SAMRecord
     */
    public static SAMRecord createArtificialSAMRecord(SAMFileHeader header, String name, int refIndex, int alignmentStart, byte[] bases, byte[] qual) {
        if (bases.length != qual.length) {
            throw new IllegalArgumentException("Passed in read string is different length then the quality array");
        }
        SAMRecord rec = createArtificialSAMRecord(header, name, refIndex, alignmentStart, bases.length);
        rec.setReadBases(Arrays.copyOf(bases, bases.length));
        rec.setBaseQualities(Arrays.copyOf(qual, qual.length));
        rec.setAttribute(SAMTag.PG.name(), new SAMReadGroupRecord("x").getId());
        if (refIndex == -1) {
            rec.setReadUnmappedFlag(true);
        }

        return rec;
    }

    /**
     * Create an artificial SAMRecord based on the parameters
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param refIndex       the reference index, i.e. what chromosome to associate it with
     * @param alignmentStart where to start the alignment
     * @param bases          the sequence of the read
     * @param qual           the qualities of the read
     * @param cigar          the cigar string of the read
     * @return the artificial SAMRecord
     */
    public static SAMRecord createArtificialSAMRecord(SAMFileHeader header, String name, int refIndex, int alignmentStart, byte[] bases, byte[] qual, String cigar) {
        SAMRecord rec = createArtificialSAMRecord(header, name, refIndex, alignmentStart, bases, qual);
        rec.setCigarString(cigar);
        return rec;
    }

    /**
     * Create an artificial SAMRecord with the following default parameters :
     * header:
     * numberOfChromosomes = 1
     * startingChromosome = 1
     * chromosomeSize = 1000000
     * read:
     * name = "default_read"
     * refIndex = 0
     * alignmentStart = 10000
     *
     * @param header SAM header for the read
     * @param bases the sequence of the read
     * @param qual  the qualities of the read
     * @param cigar the cigar string of the read
     * @return the artificial SAMRecord
     */
    public static SAMRecord createArtificialSAMRecord(final SAMFileHeader header, final byte[] bases, final byte[] qual, final String cigar) {
        return createArtificialSAMRecord(header, "default_read", 0, 10000, bases, qual, cigar);
    }

    public static SAMRecord createArtificialSAMRecord(final byte[] bases, final byte[] qual, final String cigar) {
        SAMFileHeader header = createArtificialSamHeader();
        return createArtificialSAMRecord(header, "default_read", 0, 10000, bases, qual, cigar);
    }

    public static SAMRecord createArtificialSAMRecord(final SAMFileHeader header, final Cigar cigar, final String name) {
        int length = cigar.getReadLength();
        byte base = 'A';
        byte qual = 30;
        byte [] bases = Utils.dupBytes(base, length);
        byte [] quals = Utils.dupBytes(qual, length);
        return createArtificialSAMRecord(header, name, 0, 10000, bases, quals, cigar.toString());
    }

    public static SAMRecord createArtificialSAMRecord(final SAMFileHeader header, final Cigar cigar) {
        return createArtificialSAMRecord(header, cigar, "default_read");
    }

    public static SAMRecord createArtificialSAMRecord(final Cigar cigar) {
        final SAMFileHeader header = createArtificialSamHeader();
        return createArtificialSAMRecord(header, cigar, "default_read");
    }

    /**
     * Makes a new read with a name that is unique (so that it will return false to equals(otherRead)
     */
    public static SAMRecord createUniqueArtificialSAMRecord(final Cigar cigar) {
        final SAMFileHeader header = createArtificialSamHeader();
        return createArtificialSAMRecord(header, cigar, UUID.randomUUID().toString());
    }

    public static List<GATKRead> createPair(SAMFileHeader header, String name, int readLen, int leftStart, int rightStart, boolean leftIsFirst, boolean leftIsNegative) {
        return createPair(header, name, readLen, 0, leftStart, rightStart, leftIsFirst, leftIsNegative);
    }

    public static List<GATKRead> createPair(SAMFileHeader header, String name, int readLen, int refIndex, int leftStart, int rightStart, boolean leftIsFirst, boolean leftIsNegative) {
        GATKRead left = createArtificialRead(header, name, refIndex, leftStart, readLen);
        GATKRead right = createArtificialRead(header, name, refIndex, rightStart, readLen);

        left.setIsPaired(true);
        right.setIsPaired(true);

        left.setIsProperlyPaired(true);
        right.setIsProperlyPaired(true);

        if ( leftIsFirst ) {
            left.setIsFirstOfPair();
            right.setIsSecondOfPair();
        }
        else {
            left.setIsSecondOfPair();
            right.setIsFirstOfPair();
        }

        left.setIsReverseStrand(leftIsNegative);
        left.setMateIsReverseStrand(!leftIsNegative);
        right.setIsReverseStrand(!leftIsNegative);
        right.setMateIsReverseStrand(leftIsNegative);

        left.setMatePosition(header.getSequence(refIndex).getSequenceName(), right.getStart());
        right.setMatePosition(header.getSequence(refIndex).getSequenceName(), left.getStart());

        int isize = rightStart + readLen - leftStart;
        left.setFragmentLength(isize);
        right.setFragmentLength(-isize);

        return Arrays.asList(left, right);
    }

    /**
     * Create a collection of identical artificial reads based on the parameters.  The cigar string for each
     * read will be *M, where * is the length of the read.
     *
     * Useful for testing things like positional downsampling where you care only about the position and
     * number of reads, and not the other attributes.
     *
     * @param size           number of identical reads to create
     * @param header         the SAM header to associate each read with
     * @param name           name associated with each read
     * @param refIndex       the reference index, i.e. what chromosome to associate them with
     * @param alignmentStart where to start each alignment
     * @param length         the length of each read
     *
     * @return a collection of stackSize reads all sharing the above properties
     */
    public static Collection<GATKRead> createIdenticalArtificialReads(final int size, final SAMFileHeader header, final String name, final int refIndex, final int alignmentStart, final int length ) {
        Utils.validateArg(size >= 0, "size must be non-negative");
        final Collection<GATKRead> coll = new ArrayList<>(size);
        for ( int i = 1; i <= size; i++ ) {
            coll.add(createArtificialRead(header, name, refIndex, alignmentStart, length));
        }
        return coll;
    }

    public static GATKRead createRandomRead(SAMFileHeader header, int length) {
        List<CigarElement> cigarElements = new LinkedList<>();
        cigarElements.add(new CigarElement(length, CigarOperator.M));
        Cigar cigar = new Cigar(cigarElements);
        return createArtificialRead(header, cigar);
    }

    public static GATKRead createRandomRead(int length) {
        SAMFileHeader header = createArtificialSamHeader();
        return createRandomRead(header, length);
    }

    public static GATKRead createRandomRead(SAMFileHeader header, int length, boolean allowNs) {
        byte[] quals = createRandomReadQuals(length);
        byte[] bbases = createRandomReadBases(length, allowNs);
        return createArtificialRead(bbases, quals, bbases.length + "M");
    }

    public static GATKRead createRandomRead(int length, boolean allowNs) {
        SAMFileHeader header = createArtificialSamHeader();
        return createRandomRead(header, length, allowNs);
    }

    /**
     * Create random read qualities
     *
     * @param length the length of the read
     * @return an array with randomized base qualities between 0 and 50
     */
    public static byte[] createRandomReadQuals(int length) {
        Random random = Utils.getRandomGenerator();
        byte[] quals = new byte[length];
        for (int i = 0; i < length; i++)
            quals[i] = (byte) random.nextInt(50);
        return quals;
    }

    /**
     * Create random read qualities
     *
     * @param length  the length of the read
     * @param allowNs whether or not to allow N's in the read
     * @return an array with randomized bases (A-N) with equal probability
     */
    public static byte[] createRandomReadBases(int length, boolean allowNs) {
        Random random = Utils.getRandomGenerator();
        int numberOfBases = allowNs ? 5 : 4;
        byte[] bases = new byte[length];
        for (int i = 0; i < length; i++) {
            switch (random.nextInt(numberOfBases)) {
                case 0:
                    bases[i] = 'A';
                    break;
                case 1:
                    bases[i] = 'C';
                    break;
                case 2:
                    bases[i] = 'G';
                    break;
                case 3:
                    bases[i] = 'T';
                    break;
                case 4:
                    bases[i] = 'N';
                    break;
                default:
                    throw new GATKException("Something went wrong, this is just impossible");
            }
        }
        return bases;
    }
    /**
     * create an iterator containing the specified read piles
     *
     * @param startingChr the chromosome (reference ID) to start from
     * @param endingChr   the id to end with
     * @param readCount   the number of reads per chromosome
     * @return iterator representing the specified amount of fake data
     */
    public static ArtificialReadQueryIterator mappedReadIterator(int startingChr, int endingChr, int readCount) {
        SAMFileHeader header = createArtificialSamHeader((endingChr - startingChr) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);

        return new ArtificialReadQueryIterator(startingChr, endingChr, readCount, 0, header);
    }

    /**
     * create an iterator containing the specified read piles
     *
     * @param startingChr       the chromosome (reference ID) to start from
     * @param endingChr         the id to end with
     * @param readCount         the number of reads per chromosome
     * @param unmappedReadCount the count of unmapped reads to place at the end of the iterator, like in a sorted bam file
     * @return iterator representing the specified amount of fake data
     */
    public static ArtificialReadQueryIterator mappedAndUnmappedReadIterator(int startingChr, int endingChr, int readCount, int unmappedReadCount) {
        SAMFileHeader header = createArtificialSamHeader((endingChr - startingChr) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);

        return new ArtificialReadQueryIterator(startingChr, endingChr, readCount, unmappedReadCount, header);
    }

    /**
     * Creates an artificial sam header based on the sequence dictionary dict
     *
     * @return a new sam header
     */
    public static SAMFileHeader createArtificialSamHeader(final SAMSequenceDictionary dict) {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        header.setSequenceDictionary(dict);
        return header;
    }

    /**
     * Create a pileupElement with the given insertion added to the bases.
     *
     * Assumes the insertion is prepended with one "reference" base.
     *
     * @param offsetIntoRead the offset into the read where the insertion Allele should appear.  As a reminder, the
     *                       insertion allele should have a single ref base prepend.  Must be 0 - (lengthOfRead-1)
     * @param insertionAllele the allele as you would see in a VCF for the insertion.  So, it is prepended with a ref base.  Never {@code null}
     * @param lengthOfRead the length of the artificial read.  Does not include any length differences due to the spliced indel.  Must be greater than zero.
     * @return pileupElement with an artificial read containing the insertion.
     */
    public static PileupElement createSplicedInsertionPileupElement(int offsetIntoRead, final Allele insertionAllele, final int lengthOfRead) {

        ParamUtils.isPositive(lengthOfRead, "length of read is invalid for creating an artificial read, must be greater than 0.");
        ParamUtils.inRange(offsetIntoRead, 0, lengthOfRead-1, "offset into read is invalid for creating an artificial read, must be 0-" + (lengthOfRead-1) + ".");
        Utils.nonNull(insertionAllele);

        int remainingReadLength = lengthOfRead - ((offsetIntoRead + 1) + (insertionAllele.getBases().length - 1));
        String cigarString = (offsetIntoRead + 1) + "M" + (insertionAllele.getBases().length - 1) + "I";
        if (remainingReadLength > 0) {
            cigarString += (remainingReadLength + "M");
        }

        final Cigar cigar = TextCigarCodec.decode(cigarString);
        final GATKRead gatkRead = ArtificialReadUtils.createArtificialRead(cigar);
        final PileupElement pileupElement = PileupElement.createPileupForReadAndOffset(gatkRead, offsetIntoRead);

        // Splice in that insertion.
        final byte[] bases = gatkRead.getBases();
        final int newReadLength = lengthOfRead + insertionAllele.getBases().length - 1;
        final byte[] destBases = new byte[newReadLength];
        final byte[] basesToInsert = ArrayUtils.subarray(insertionAllele.getBases(), 1, insertionAllele.getBases().length);
        System.arraycopy(bases, 0, destBases, 0, offsetIntoRead);

        // Make sure that the one prepended "reference base" matches the input.
        destBases[offsetIntoRead] = insertionAllele.getBases()[0];

        System.arraycopy(basesToInsert, 0, destBases, offsetIntoRead+1, basesToInsert.length);

        if ((offsetIntoRead + 1) < lengthOfRead) {
            System.arraycopy(bases, offsetIntoRead + 1, destBases, offsetIntoRead + basesToInsert.length + 1, bases.length - 1 - offsetIntoRead);
        }

        gatkRead.setBases(destBases);

        return pileupElement;
    }

    /** See {@link ArtificialReadUtils#createSplicedInsertionPileupElement}, except that this method returns a
     * pileup element containing the specified deletion.
     *
     * @param offsetIntoRead See {@link ArtificialReadUtils#createSplicedInsertionPileupElement}
     * @param referenceAllele  the reference allele as you would see in a VCF for the deletion.
     *                         In other words, it is the deletion prepended with a single ref base.  Never {@code null}
     * @param lengthOfRead See {@link ArtificialReadUtils#createSplicedInsertionPileupElement}
     * @return pileupElement with an artificial read containing the deletion.
     */
    public static PileupElement createSplicedDeletionPileupElement(int offsetIntoRead, final Allele referenceAllele, final int lengthOfRead) {
        ParamUtils.isPositive(lengthOfRead, "length of read is invalid for creating an artificial read, must be greater than 0.");
        ParamUtils.inRange(offsetIntoRead, 0, lengthOfRead-1, "offset into read is invalid for creating an artificial read, must be 0-" + (lengthOfRead-1) + ".");
        Utils.nonNull(referenceAllele);

        // Do not include the prepended "ref"
        final int numberOfSpecifiedBasesToDelete = referenceAllele.getBases().length - 1;
        final int numberOfBasesToActuallyDelete = Math.min(numberOfSpecifiedBasesToDelete, lengthOfRead - offsetIntoRead - 1);

        final int newReadLength = lengthOfRead - numberOfBasesToActuallyDelete;

        String cigarString = (offsetIntoRead + 1) + "M";

        if (numberOfBasesToActuallyDelete > 0) {
            cigarString += numberOfBasesToActuallyDelete + "D";
        }
        final int remainingBases = lengthOfRead - (offsetIntoRead + 1) - numberOfBasesToActuallyDelete;
        if (remainingBases > 0) {
            cigarString += remainingBases + "M";
        }

        final Cigar cigar = TextCigarCodec.decode(cigarString);
        final GATKRead gatkRead = ArtificialReadUtils.createArtificialRead(cigar);
        final PileupElement pileupElement = PileupElement.createPileupForReadAndOffset(gatkRead, offsetIntoRead);

        // The Cigar string has basically already told the initial generation of a read to delete bases.
        final byte[] bases = gatkRead.getBases();

        // Make sure that the one prepended "reference base" matches the input.
        bases[offsetIntoRead] = referenceAllele.getBases()[0];

        gatkRead.setBases(bases);

        return pileupElement;
    }

    /**
     * See {@link ArtificialReadUtils#createSplicedInsertionPileupElement}, except that this method returns a
     * pileup element containing base-by-base replacement.  As a result, the length of the read will not change.
     *
     * @param offsetIntoRead See {@link ArtificialReadUtils#createSplicedInsertionPileupElement}
     * @param newAllele The new bases that should be in the read at the specified position.  If this allele causes the
     *                  replacement to extend beyond the end of the read
     *                  (i.e. offsetIntoRead + length(newAllele) is greater than length of read),
     *                  the replacement will be truncated.
     * @param lengthOfRead See {@link ArtificialReadUtils#createSplicedInsertionPileupElement}
     * @return pileupElement with an artificial read containing the new bases specified by te given allele.
     */
    public static PileupElement createNonIndelPileupElement(final int offsetIntoRead, final Allele newAllele, final int lengthOfRead) {
        ParamUtils.isPositive(lengthOfRead, "length of read is invalid for creating an artificial read, must be greater than 0.");
        ParamUtils.inRange(offsetIntoRead, 0, lengthOfRead-1, "offset into read is invalid for creating an artificial read, must be 0-" + (lengthOfRead-1) + ".");
        Utils.nonNull(newAllele);

        final String cigarString = lengthOfRead + "M";

        final Cigar cigar = TextCigarCodec.decode(cigarString);
        final GATKRead gatkRead = ArtificialReadUtils.createArtificialRead(cigar);
        final byte[] newBases = gatkRead.getBases();
        final int upperBound = Math.min(offsetIntoRead + newAllele.getBases().length, lengthOfRead);
        for (int i = offsetIntoRead; i < upperBound; i++) {
            newBases[i] = newAllele.getBases()[i - offsetIntoRead];
        }
        gatkRead.setBases(newBases);
        return PileupElement.createPileupForReadAndOffset(gatkRead, offsetIntoRead);
    }
}
