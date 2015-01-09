package org.broadinstitute.hellbender.utils.sam;

/**
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

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.GATKSAMIterator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class ArtificialSAMUtils {
    public static final int DEFAULT_READ_LENGTH = 50;

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
        header.setSortOrder(htsjdk.samtools.SAMFileHeader.SortOrder.coordinate);
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
     * Creates an artificial sam header with standard test parameters
     *
     * @return the sam header
     */
    public static SAMFileHeader createArtificialSamHeader() {
        return createArtificialSamHeader(1, 1, 1000000);
    }

    /**
     * Create an artificial read based on the parameters.  The cigar string will be *M, where * is the length of the read
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param refIndex       the reference index, i.e. what chromosome to associate it with
     * @param alignmentStart where to start the alignment
     * @param length         the length of the read
     * @return the artificial read
     */
    public static SAMRecord createArtificialRead(SAMFileHeader header, String name, int refIndex, int alignmentStart, int length) {
        if ((refIndex == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && alignmentStart != SAMRecord.NO_ALIGNMENT_START) ||
                (refIndex != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && alignmentStart == SAMRecord.NO_ALIGNMENT_START))
            throw new IllegalArgumentException("Invalid alignment start for artificial read, start = " + alignmentStart);
        SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceIndex(refIndex);
        record.setAlignmentStart(alignmentStart);
        List<CigarElement> elements = new ArrayList<CigarElement>();
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
     * Create an artificial read based on the parameters.  The cigar string will be *M, where * is the length of the read
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param refIndex       the reference index, i.e. what chromosome to associate it with
     * @param alignmentStart where to start the alignment
     * @param bases          the sequence of the read
     * @param qual           the qualities of the read
     * @return the artificial read
     */
    public static SAMRecord createArtificialRead(SAMFileHeader header, String name, int refIndex, int alignmentStart, byte[] bases, byte[] qual) {
        if (bases.length != qual.length) {
            throw new IllegalArgumentException("Passed in read string is different length then the quality array");
        }
        SAMRecord rec = createArtificialRead(header, name, refIndex, alignmentStart, bases.length);
        rec.setReadBases(bases);
        rec.setBaseQualities(qual);
        rec.setAttribute(SAMTag.PG.name(), new SAMReadGroupRecord("x").getId());
        if (refIndex == -1) {
            rec.setReadUnmappedFlag(true);
        }

        return rec;
    }

    /**
     * Create an artificial read based on the parameters
     *
     * @param header         the SAM header to associate the read with
     * @param name           the name of the read
     * @param refIndex       the reference index, i.e. what chromosome to associate it with
     * @param alignmentStart where to start the alignment
     * @param bases          the sequence of the read
     * @param qual           the qualities of the read
     * @param cigar          the cigar string of the read
     * @return the artificial read
     */
    public static SAMRecord createArtificialRead(SAMFileHeader header, String name, int refIndex, int alignmentStart, byte[] bases, byte[] qual, String cigar) {
        SAMRecord rec = createArtificialRead(header, name, refIndex, alignmentStart, bases, qual);
        rec.setCigarString(cigar);
        return rec;
    }

    /**
     * Create an artificial read with the following default parameters :
     * header:
     * numberOfChromosomes = 1
     * startingChromosome = 1
     * chromosomeSize = 1000000
     * read:
     * name = "default_read"
     * refIndex = 0
     * alignmentStart = 10000
     *
     * @param bases the sequence of the read
     * @param qual  the qualities of the read
     * @param cigar the cigar string of the read
     * @return the artificial read
     */
    public static SAMRecord createArtificialRead(byte[] bases, byte[] qual, String cigar) {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader();
        return ArtificialSAMUtils.createArtificialRead(header, "default_read", 0, 10000, bases, qual, cigar);
    }

    public static SAMRecord createArtificialRead(Cigar cigar) {
        int length = cigar.getReadLength();
        byte [] base = {'A'};
        byte [] qual = {30};
        byte [] bases = Utils.arrayFromArrayWithLength(base, length);
        byte [] quals = Utils.arrayFromArrayWithLength(qual, length);
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader();
        return ArtificialSAMUtils.createArtificialRead(header, "default_read", 0, 10000, bases, quals, cigar.toString());
    }


    public final static List<SAMRecord> createPair(SAMFileHeader header, String name, int readLen, int leftStart, int rightStart, boolean leftIsFirst, boolean leftIsNegative) {
        SAMRecord left = ArtificialSAMUtils.createArtificialRead(header, name, 0, leftStart, readLen);
        SAMRecord right = ArtificialSAMUtils.createArtificialRead(header, name, 0, rightStart, readLen);

        left.setReadPairedFlag(true);
        right.setReadPairedFlag(true);

        left.setProperPairFlag(true);
        right.setProperPairFlag(true);

        left.setFirstOfPairFlag(leftIsFirst);
        right.setFirstOfPairFlag(!leftIsFirst);

        left.setReadNegativeStrandFlag(leftIsNegative);
        left.setMateNegativeStrandFlag(!leftIsNegative);
        right.setReadNegativeStrandFlag(!leftIsNegative);
        right.setMateNegativeStrandFlag(leftIsNegative);

        left.setMateAlignmentStart(right.getAlignmentStart());
        right.setMateAlignmentStart(left.getAlignmentStart());

        left.setMateReferenceIndex(0);
        right.setMateReferenceIndex(0);

        int isize = rightStart + readLen - leftStart;
        left.setInferredInsertSize(isize);
        right.setInferredInsertSize(-isize);

        return Arrays.asList(left, right);
    }

    public static SAMRecord createRandomRead(int length) {
        List<CigarElement> cigarElements = new LinkedList<>();
        cigarElements.add(new CigarElement(length, CigarOperator.M));
        Cigar cigar = new Cigar(cigarElements);
        return ArtificialSAMUtils.createArtificialRead(cigar);
    }

    /**
     * create an iterator containing the specified read piles
     *
     * @param startingChr the chromosome (reference ID) to start from
     * @param endingChr   the id to end with
     * @param readCount   the number of reads per chromosome
     * @return GATKSAMIterator representing the specified amount of fake data
     */
    public static GATKSAMIterator mappedReadIterator(int startingChr, int endingChr, int readCount) {
        SAMFileHeader header = createArtificialSamHeader((endingChr - startingChr) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);

        return new ArtificialSAMQueryIterator(startingChr, endingChr, readCount, 0, header);
    }

    /**
     * create an iterator containing the specified read piles
     *
     * @param startingChr       the chromosome (reference ID) to start from
     * @param endingChr         the id to end with
     * @param readCount         the number of reads per chromosome
     * @param unmappedReadCount the count of unmapped reads to place at the end of the iterator, like in a sorted bam file
     * @return GATKSAMIterator representing the specified amount of fake data
     */
    public static GATKSAMIterator mappedAndUnmappedReadIterator(int startingChr, int endingChr, int readCount, int unmappedReadCount) {
        SAMFileHeader header = createArtificialSamHeader((endingChr - startingChr) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);

        return new ArtificialSAMQueryIterator(startingChr, endingChr, readCount, unmappedReadCount, header);
    }

    /**
     * Creates an artificial sam header based on the sequence dictionary dict
     *
     * @return a new sam header
     */
    public static SAMFileHeader createArtificialSamHeader(final SAMSequenceDictionary dict) {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(htsjdk.samtools.SAMFileHeader.SortOrder.coordinate);
        header.setSequenceDictionary(dict);
        return header;
    }
}
