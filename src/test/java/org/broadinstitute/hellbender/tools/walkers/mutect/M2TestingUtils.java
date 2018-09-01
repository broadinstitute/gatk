package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class M2TestingUtils {
    private static final String DEFAULT_READ_GROUP_NAME = "READ_GROUP_1";
    /**
     *
     * Reference at this position looks like this:
     * coord.:    99,990             100,000             100,010             100,020
     *              |                   |                   |                   |
     * reference:   T C A T C A C A C T C A C T A A G C A C A C A G A G A A T A A T
     * alt read:    T C A T C A C A C T C T C T G A C C A A A G A G A A T A A T A A
     */
    private static final int DEFAULT_ALIGNMENT_START = 99_991;
    private static final int DEFAULT_CHROM_INDEX = 0;
    private static final byte DEFAULT_BASEQ = 26;
    private static final int DEFAULT_MAPQ = 60;
    private static final int DEFAULT_NUM_CHROMOSOMES = 1;
    private static final int DEFAULT_STARTING_CHROM = 1;
    private static final int DEFAULT_CHROMOSOME_SIZE = 1_000_000;
    public static final byte[] DEFAULT_REF_BASES = "CATCACACTCACTAAGCACACAGAGAATAAT".getBytes();
    // Bases for C->T SNP at position 100,000                *
    public static final byte[] DEFAULT_ALT_BASES = "CATCACACTTACTAAGCACACAGAGAATAAT".getBytes();
    public static final int DEFAULT_SNP_POSITION = 100_000;
    public static final String DEFAULT_SAMPLE_NAME = "sample1";

    public static SAMFileGATKReadWriter getBareBonesSamWriter(final File samFile, final SAMFileHeader samHeader) {
        return new SAMFileGATKReadWriter(ReadUtils.createCommonSAMWriter(samFile, null, samHeader, true, true, false));
    }

    // TODO: might have to address read name collisions
    public static List<GATKRead> createReads(final int numReads, final byte[] bases, final SAMFileHeader samHeader,
                                             final byte baseQuality){
        final int readLength = bases.length;
        final byte[] quals = new byte[readLength];
        Arrays.fill(quals, baseQuality);

        final List<GATKRead> reads = new ArrayList<>(numReads);
        for (int i = 0; i < numReads; i++) {
            final String cigar = bases.length + "M";
            final GATKRead read = ArtificialReadUtils.createArtificialRead(samHeader, "Read" + i, DEFAULT_CHROM_INDEX,
                    DEFAULT_ALIGNMENT_START, bases, quals, cigar);
            read.setReadGroup(DEFAULT_READ_GROUP_NAME);
            read.setMappingQuality(DEFAULT_MAPQ);
            read.setIsFirstOfPair();
            read.setIsReverseStrand(i % 2 == 0);
            reads.add(read);
        }

        return reads;
    }

    public static List<GATKRead> createReads(final int numReads, final byte[] bases, final SAMFileHeader samHeader){
        return createReads(numReads, bases, samHeader, DEFAULT_BASEQ);
    }

    public static SAMFileHeader createSamHeader(){
        final SAMFileHeader samHeader = ArtificialReadUtils.createArtificialSamHeader(
                DEFAULT_NUM_CHROMOSOMES, DEFAULT_STARTING_CHROM, DEFAULT_CHROMOSOME_SIZE);
        final SAMReadGroupRecord readGroupRecord = new SAMReadGroupRecord(DEFAULT_READ_GROUP_NAME);
        readGroupRecord.setSample(DEFAULT_SAMPLE_NAME);
        samHeader.addReadGroup(readGroupRecord);
        return samHeader;
    }

}
