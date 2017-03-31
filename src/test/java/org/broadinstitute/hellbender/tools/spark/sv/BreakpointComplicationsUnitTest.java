package org.broadinstitute.hellbender.tools.spark.sv;


import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

public class BreakpointComplicationsUnitTest {

    @Test
    public void testGetHomology() {
        final byte[] contigSequence = "ATCGATCGAAAAGCTAGCTA".getBytes();

        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 1, 12), TextCigarCodec.decode("12M8S"), true, 60, 1, 1, 12);            // dummy test data, almost guaranteed to be non-factual
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 101, 112), TextCigarCodec.decode("8H12M"), false, 60, 1, 9, 20);    // dummy test data, almost guaranteed to be non-factual

        Assert.assertEquals(BreakpointComplications.getHomology(region1, region2, contigSequence), "AAAA");

        final AlignmentRegion region3 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 1, 12), TextCigarCodec.decode("8M"), true, 60, 1, 1, 8);            // dummy test data, almost guaranteed to be non-factual
        final AlignmentRegion region4 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 101, 112), TextCigarCodec.decode("8M"), false, 60, 1, 13, 20);    // dummy test data, almost guaranteed to be non-factual

        Assert.assertTrue(BreakpointComplications.getHomology(region3, region4, contigSequence).isEmpty());
    }

    @Test
    public void testGetInsertedSequence() {
        final byte[] contigSequence = "GACGAACGATTTGACTTTAATATGAAATGTTTTATGTGGGCTATAAAATTATCCAAACTCGACACAGGACATTTTGAGCTTATTTCCAAATCATCTGGCCTTCATCTACCCACTGGAACTATTACTCTGCTGGGTCCTCATGGAAACATATCTTTCAGCCCTAACAATGAGACTACAGACATCTACGTCCCCAACACAACAGCTAAAAAGCAGTAGAATGTCAGAAAGGCTATCCACTTAGCCCTTGGCTGACAGGCCCCACTGAGCATCCTTTGCGAAGTCCATTTACTAGCTAATTCATAATTTACACAAGGCATTCAGACATAGCAGCTAAGATATAAAACATTTATCAACACAGGGACTAGTTTGTCATTTTAAAATAATTATGTTTAAGTAAGCCAATAAAGTCTATCTTCTCCAATTTACTTATTGAGCTTTATGAGGCAATTTAAGTCCCGATTTTGGGGGGTATGTATGAAAGGAGAGCATGGAAATGCCATTTGCTCCCTGAAGTTTTTATCTTTTTTTTTTTGAGATAGAGTCTTGTGTTTTCTGTGGAGTACATGAGTATGCATCAAAGCTAACAACGCCCACTGCCCTGTTAGTCAAATACCTTTGA".getBytes();
        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("8", 118873207, 118873739), TextCigarCodec.decode("532M87S"), true, 60, 0, 1, 532);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 175705642, 175705671), TextCigarCodec.decode("518S29M72S"), false, 3, 0, 519, 547);
        final AlignmentRegion region3 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 118875262, 118875338), TextCigarCodec.decode("543S76M"), false, 60, 0, 544, 619);

        Assert.assertTrue(BreakpointComplications.getInsertedSequence(region3, region1, contigSequence).isEmpty());
        Assert.assertEquals(BreakpointComplications.getInsertedSequence(region1, region3, contigSequence), "GAGATAGAGTC");

        Assert.assertTrue(BreakpointComplications.getInsertedSequence(region2, region1, contigSequence).isEmpty() && BreakpointComplications.getInsertedSequence(region1, region2, contigSequence).isEmpty());
    }

    // -----------------------------------------------------------------------------------------------
    // Tests for CIGAR extraction on tandem duplications
    // -----------------------------------------------------------------------------------------------
    @Test
    public void testExtractCigar() {

        final int contigTotalLength = 355;

        // forward strand
        final AlignmentRegion region1 = new AlignmentRegion(AlignmentRegion.DUMMY_ASM_ID, "1", new SimpleInterval("1", 1000001, 1000125),
                TextCigarCodec.decode("5H10S15M20D25M30D35M260S5H"),
                true, 60, 0, 16, 75);
        final AlignmentRegion region2 = new AlignmentRegion(AlignmentRegion.DUMMY_ASM_ID, "1", new SimpleInterval("1", 1000041, 1000145),
                TextCigarCodec.decode("5H185S45M30I55M20I5M10S5H"),
                true, 60, 0, 191, 340);

        final Cigar cigar1 = BreakpointComplications.extractCigarForTandup(region1, 1000125, 1000041);
        Assert.assertEquals(cigar1, TextCigarCodec.decode("20M30D35M"));
        final Cigar cigar2 = BreakpointComplications.extractCigarForTandup(region2, 1000125, 1000041);
        Assert.assertEquals(cigar2, TextCigarCodec.decode("45M30I40M"));

        // reverse strand
        final AlignmentRegion region3 = new AlignmentRegion(AlignmentRegion.DUMMY_ASM_ID, "1", region2.referenceInterval,
                CigarUtils.invertCigar(region2.cigarAlong5to3DirectionOfContig),
                false, 60, 0, contigTotalLength-region2.endInAssembledContig+1, contigTotalLength-region2.startInAssembledContig+1);
        final AlignmentRegion region4 = new AlignmentRegion(AlignmentRegion.DUMMY_ASM_ID, "1", region1.referenceInterval,
                CigarUtils.invertCigar(region1.cigarAlong5to3DirectionOfContig),
                false, 60, 0, contigTotalLength-region1.endInAssembledContig+1, contigTotalLength-region1.startInAssembledContig+1);

        final Cigar cigar3 = BreakpointComplications.extractCigarForTandup(region3, 1000125, 1000041);
        Assert.assertEquals(CigarUtils.invertCigar(cigar3), cigar2);
        final Cigar cigar4 = BreakpointComplications.extractCigarForTandup(region4, 1000125, 1000041);
        Assert.assertEquals(CigarUtils.invertCigar(cigar4), cigar1);
    }
}
