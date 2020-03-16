package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.engine.filters.UMIReadFilter;
import org.broadinstitute.hellbender.utils.read.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class M2TestingUtils {
    public static final String DEFAULT_READ_GROUP_NAME = "READ_GROUP_1";
    /**
     *
     * hg19 reference at this position in chr1 looks like this:
     * coord.:    99,990             100,000             100,010             100,020
     *              |                   |                   |                   |
     * reference:   T C A T C A C A C T C A C T A A G C A C A C A G A G A A T A A T G T C
     * alt read:    T C A T C A C A C T C T C T G A C C A A A G A G A A T A A T A A G T C
     *
     * And its mate 100 bps upstream
     coord.:     100,090             100,100             100,110             100,120
     *              |                   |                   |                   |
     * reference:   A G C A T C G C C A T G C C T A G T A C A G A C T C T C C C T
     *
     */
    private static final int DEFAULT_ALIGNMENT_START = 99_991;
    private static final int MATE_START = 100_090;
    private static final int DEFAULT_CHROM_INDEX = 0;
    private static final int DEFAULT_MAPQ = 60;
    private static final int DEFAULT_NUM_CHROMOSOMES = 1;
    private static final int DEFAULT_STARTING_CHROM = 1;
    private static final int DEFAULT_CHROMOSOME_SIZE = 1_000_000;

    // index:   1234567890123456789012345678901
    // ref:     CATCACACTCACTAAGCACACAGAGAATAAT
    // AC del:  CATCACACTCACTAAGC--ACAGAGAATAATGT    cigar = "17M2D13M"
    // ACA del: CATCACACTCACTAAGC---CAGAGAATAATGTC   cigar = "17M3D13M"


    // index:   1234567890123456789012345678901
    // ref:     CATCACACTCACTAAGC---ACACAGAGAATAAT
    // TTT ins: CATCACACTCACTAAGCTTTACACAGAGAAT         cigar = "17M3I10M"

    public static final byte[] DEFAULT_REF_BASES = "CATCACACTCACTAAGCACACAGAGAATAAT".getBytes();
    public static final byte[] MATE_REF_BASES =    "AGCATCGCCATGCCTAGTACAGACTCTCCCT".getBytes();
    public static final byte[] AC_DELETION =       "CATCACACTCACTAAGCACAGAGAATAATGT".getBytes();
    public static final byte[] ACA_DELETION =      "CATCACACTCACTAAGCCAGAGAATAATGTC".getBytes();
    public static final byte[] TTT_INSERTION =     "CATCACACTCACTAAGCTTTACACAGAGAAT".getBytes();
    public static final int DEFAULT_REF_READ_LENGTH = DEFAULT_REF_BASES.length;
    public static final String CIGAR_REF = "31M";
    public static final String CIGAR_AC_DELETION = "17M2D14M";
    public static final String CIGAR_ACA_DELETION= "17M3D14M";
    public static final String CIGAR_TTT_INSERTION = "17M3I11M";

    public static final String DEFAULT_UMI = "ACT-GGG";
    public static final String DEFAULT_MATE_UMI = "GGG-ACT";

    public static final int FRAGMENT_LENGTH = MATE_START + DEFAULT_REF_READ_LENGTH - DEFAULT_ALIGNMENT_START;


    public static SAMFileGATKReadWriter getBareBonesSamWriter(final File samFile, final SAMFileHeader samHeader) {
        return new SAMFileGATKReadWriter(ReadUtils.createCommonSAMWriter(samFile, null, samHeader, true, true, false));
    }

    public static List<GATKRead> createReads(final int numReads, final byte[] bases, final SAMFileHeader samHeader,
                                             final byte baseQuality, final String namePrefix) {
        return createReads(numReads, bases, samHeader, baseQuality, namePrefix, bases.length + "M", "AAA-AAA", false, true);
    }

    public static List<GATKRead> createReads(final int numReads, final byte[] bases, final SAMFileHeader samHeader,
                                             final byte baseQuality, final String namePrefix, final String cigar,
                                             final String umi, final boolean withMate, final boolean f1r2){
        final int readLength = bases.length;

        final List<GATKRead> reads = new ArrayList<>(numReads);
        final List<GATKRead> mates = new ArrayList<>(numReads);
        for (int i = 0; i < numReads; i++) {
            // final String cigar = CigarUtils.calculateCigar(DEFAULT_REF_BASES, , new SmithWatermanAligner(),
            final byte[] quals = new byte[readLength];
            Arrays.fill(quals, baseQuality);
            final GATKRead read = ArtificialReadUtils.createArtificialRead(samHeader, namePrefix + i + f1r2, DEFAULT_CHROM_INDEX,
                    DEFAULT_ALIGNMENT_START, bases, quals, cigar);
            read.setReadGroup(DEFAULT_READ_GROUP_NAME);
            read.setMappingQuality(DEFAULT_MAPQ);
            read.setAttribute(UMIReadFilter.UMI_TAG, M2TestingUtils.DEFAULT_UMI);
            if (f1r2) { read.setIsFirstOfPair(); } else { read.setIsSecondOfPair(); }

            if (withMate){
                read.setMatePosition(Integer.toString(DEFAULT_STARTING_CHROM), MATE_START);
                read.isPaired();
                read.setMateIsReverseStrand(true);
                read.setIsProperlyPaired(true);
                read.setAttribute("MQ", 60); // mate quality!
                read.setAttribute("MC", CIGAR_REF);

                final GATKRead mate = ArtificialReadUtils.createArtificialRead(samHeader, namePrefix + i + f1r2, DEFAULT_CHROM_INDEX,
                        MATE_START, MATE_REF_BASES, quals, CIGAR_REF);
                mate.setReadGroup(DEFAULT_READ_GROUP_NAME);
                mate.setMappingQuality(DEFAULT_MAPQ);
                if (f1r2) { mate.setIsSecondOfPair(); } else { mate.setIsFirstOfPair(); }
                mate.setMatePosition(Integer.toString(DEFAULT_STARTING_CHROM), DEFAULT_ALIGNMENT_START);
                mate.setAttribute(UMIReadFilter.UMI_TAG, M2TestingUtils.DEFAULT_MATE_UMI);
                mate.isPaired();
                mate.setIsReverseStrand(true);
                mate.setIsProperlyPaired(true);
                mate.setAttribute("MQ", 60);
                mate.setAttribute("MC", cigar);

                read.setFragmentLength(FRAGMENT_LENGTH);
                mate.setFragmentLength(-FRAGMENT_LENGTH); // See sam spec TLEN for how insert size is encoded
                mates.add(mate);
            } else {
                read.setMateIsUnplaced();
            }

            reads.add(read);
        }

        reads.addAll(mates);
        return reads;
    }

    public static SAMFileHeader createSamHeader(final String sampleName){
        final SAMFileHeader samHeader = ArtificialReadUtils.createArtificialSamHeader(
                DEFAULT_NUM_CHROMOSOMES, DEFAULT_STARTING_CHROM, DEFAULT_CHROMOSOME_SIZE);
        final SAMReadGroupRecord readGroupRecord = new SAMReadGroupRecord(DEFAULT_READ_GROUP_NAME);
        readGroupRecord.setSample(sampleName);
        samHeader.addReadGroup(readGroupRecord);
        return samHeader;
    }

}
