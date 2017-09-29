package org.broadinstitute.hellbender.tools.spark.sv.discovery;


import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

public class BreakpointComplicationsUnitTest {

    @Test(groups = "sv")
    public void testGetHomology() {
        final byte[] contigSequence = "ATCGATCGAAAAGCTAGCTA".getBytes();

        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("1", 1, 12), 1, 12, TextCigarCodec.decode("12M8S"), true, 60, 1, 100, false, false);            // dummy test data, almost guaranteed to be non-factual
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("1", 101, 112), 9, 20, TextCigarCodec.decode("8H12M"), false, 60, 1, 100, false, false);    // dummy test data, almost guaranteed to be non-factual

        Assert.assertEquals(BreakpointComplications.getHomology(region1, region2, contigSequence), "AAAA");

        final AlignmentInterval region3 = new AlignmentInterval(new SimpleInterval("1", 1, 12), 1, 8, TextCigarCodec.decode("8M"), true, 60, 1, 100, false, false);            // dummy test data, almost guaranteed to be non-factual
        final AlignmentInterval region4 = new AlignmentInterval(new SimpleInterval("1", 101, 112), 13, 20, TextCigarCodec.decode("8M"), false, 60, 1, 100, false, false);    // dummy test data, almost guaranteed to be non-factual

        Assert.assertTrue(BreakpointComplications.getHomology(region3, region4, contigSequence).isEmpty());
    }

    @Test(groups = "sv")
    public void testGetInsertedSequence() {
        final byte[] contigSequence = "GACGAACGATTTGACTTTAATATGAAATGTTTTATGTGGGCTATAAAATTATCCAAACTCGACACAGGACATTTTGAGCTTATTTCCAAATCATCTGGCCTTCATCTACCCACTGGAACTATTACTCTGCTGGGTCCTCATGGAAACATATCTTTCAGCCCTAACAATGAGACTACAGACATCTACGTCCCCAACACAACAGCTAAAAAGCAGTAGAATGTCAGAAAGGCTATCCACTTAGCCCTTGGCTGACAGGCCCCACTGAGCATCCTTTGCGAAGTCCATTTACTAGCTAATTCATAATTTACACAAGGCATTCAGACATAGCAGCTAAGATATAAAACATTTATCAACACAGGGACTAGTTTGTCATTTTAAAATAATTATGTTTAAGTAAGCCAATAAAGTCTATCTTCTCCAATTTACTTATTGAGCTTTATGAGGCAATTTAAGTCCCGATTTTGGGGGGTATGTATGAAAGGAGAGCATGGAAATGCCATTTGCTCCCTGAAGTTTTTATCTTTTTTTTTTTGAGATAGAGTCTTGTGTTTTCTGTGGAGTACATGAGTATGCATCAAAGCTAACAACGCCCACTGCCCTGTTAGTCAAATACCTTTGA".getBytes();
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("8", 118873207, 118873739), 1, 532, TextCigarCodec.decode("532M87S"), true, 60, 0, 100, false, false);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("1", 175705642, 175705671), 519, 547, TextCigarCodec.decode("518S29M72S"), false, 3, 0, 100, false, false);
        final AlignmentInterval region3 = new AlignmentInterval(new SimpleInterval("1", 118875262, 118875338), 544, 619, TextCigarCodec.decode("543S76M"), false, 60, 0, 100, false, false);

        Assert.assertTrue(BreakpointComplications.getInsertedSequence(region3, region1, contigSequence).isEmpty());
        Assert.assertEquals(BreakpointComplications.getInsertedSequence(region1, region3, contigSequence), "GAGATAGAGTC");

        Assert.assertTrue(BreakpointComplications.getInsertedSequence(region2, region1, contigSequence).isEmpty() && BreakpointComplications.getInsertedSequence(region1, region2, contigSequence).isEmpty());
    }

    @Test(groups = "sv")
    public void testIsLikelyInvertedDuplication() {

        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(21, 1, 46709983);

        final SAMRecord one =
                ArtificialReadUtils.createArtificialRead(header, "asm024381:tig00001", 20, 6070057, "TGAAATGTATGTGTGAGGTATGCAGTATGTGTGTGAGGTAGTGTGCGATGTGTGTGTAGTGTATGTGGTTGTGTGAGGTATGTGGGGTGTGAGGAATGTATTGTGTATGTGTGATATATACTGATTGTGTGTAAGGGATGTGGAGTGTGTGGCATGTGTGTAAGGTAGGTGTGTGTTGTGTATATGTGAGCTGTATAGTGTCGGGGGGGTGTGAGGTATGTGGTGTATGTTATGTTTGAGATCTAGTGTGTGTGTATGGTGTGTGTGGGAGGTATGTGGGGTGTGTGGTGTGTGGTGTGTATGAGGTATGTAGTGTGAGGTGTGTGATGTGTAGTGTGTGGTGTGGGGTATGTGGTGTATGTGTGAAGTATGTGTTGTGTGATGTGTGGGTGATATTTGGTGCCGTGTGTGTGGTATATGGTGTGTGGTATGAGGTGTGTAGTGTGATATGTGTGGTGTGTAATATGTGGTGTGTGTGTGTGTGTGATATATGGTGTGTGTGGTGTTATGATGTGTGTTGTGAGGTATGTGGTGTCTGTGTGTGATATGTGATTTGGGTGTGAGGTGTGTGTGGTGTGGCGTGTGGTGTGTGTGATGTGATGTGTGTGTGACATGGGGTGGTGCGTGGTGTGGTGTGTGTGGTATGTGGTGGTTGGTGTGTATGTGGTGAGTGAGGGGTGTGTGGTGTGGGTGGTGTGTGTGGTGTGTGTGGTTTGTGGTGTGTGTGGTTTGTGGTGTGTGGTATGTGGTGTGTTGTGTGTGGTTTGTGGTATGGTGTGTGTGGTATGGTTGTGTGTGGTGTGGTGTGTGCTGTGTGTATGGTTTGTGGTGTGTGTGGTGTGT".getBytes(),
                        ArtificialReadUtils.createRandomReadQuals(843), "502M341S").convertToSAMRecord(header);
        one.setMappingQuality(60);
        one.setAttribute("NM", 0);
        one.setAttribute("AS", 502);

        final SAMRecord two =
                ArtificialReadUtils.createArtificialRead(header, "asm024381:tig00001", 20, 43467994, "ACACACCACACACACCACAAACCATACACACAGCACACACCACACCACACACAACCATACCACACACACCATACCACAAACCACACACAACACACCACATACCACACACCACAAACCACACACACCACAAACCACACACACCACACACACCACCCACACCACACACC".getBytes(),
                        ArtificialReadUtils.createRandomReadQuals(167), "167M676H").convertToSAMRecord(header);
        two.setSupplementaryAlignmentFlag(true);
        two.setMappingQuality(60);
        two.setReadNegativeStrandFlag(true);
        two.setAttribute("NM", 0);
        two.setAttribute("AS", 167);

        final AlignmentInterval intervalOne = new AlignmentInterval(one);
        final AlignmentInterval intervalTwo = new AlignmentInterval(two);

        final AlignedContig contig = new AlignedContig("asm024381:tig00001", one.getReadBases(),
                Arrays.asList(intervalOne, intervalTwo), false);

        Assert.assertFalse( BreakpointComplications.isLikelyInvertedDuplication(intervalOne, intervalTwo) );
    }

    // -----------------------------------------------------------------------------------------------
    // Tests for CIGAR extraction on tandem duplications
    // -----------------------------------------------------------------------------------------------
    @Test(groups = "sv")
    public void testExtractCigarForSimleTandup() {

        final int contigTotalLength = 355;

        // forward strand
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("1", 1000001, 1000125), 16, 75,
                TextCigarCodec.decode("5H10S15M20D25M30D35M260S5H"),
                true, 60, 0, 100, false, false);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("1", 1000041, 1000145), 191, 340,
                TextCigarCodec.decode("5H185S45M30I55M20I5M10S5H"),
                true, 60, 0, 100, false, false);


        final Cigar cigar1 = BreakpointComplications.extractCigarForTandup(region1, 1000125, 1000041);
        Assert.assertEquals(cigar1, TextCigarCodec.decode("20M30D35M"));
        final Cigar cigar2 = BreakpointComplications.extractCigarForTandup(region2, 1000125, 1000041);
        Assert.assertEquals(cigar2, TextCigarCodec.decode("45M30I40M"));

        // reverse strand
        final AlignmentInterval region3 = new AlignmentInterval(region2.referenceSpan, contigTotalLength-region2.endInAssembledContig+1, contigTotalLength-region2.startInAssembledContig+1,
                CigarUtils.invertCigar(region2.cigarAlong5to3DirectionOfContig),
                false, 60, 0, 100, false, false);
        final AlignmentInterval region4 = new AlignmentInterval(region1.referenceSpan, contigTotalLength-region1.endInAssembledContig+1, contigTotalLength-region1.startInAssembledContig+1,
                CigarUtils.invertCigar(region1.cigarAlong5to3DirectionOfContig),
                false, 60, 0, 100, false, false);

        final Cigar cigar3 = BreakpointComplications.extractCigarForTandup(region3, 1000125, 1000041);
        Assert.assertEquals(CigarUtils.invertCigar(cigar3), cigar2);
        final Cigar cigar4 = BreakpointComplications.extractCigarForTandup(region4, 1000125, 1000041);
        Assert.assertEquals(CigarUtils.invertCigar(cigar4), cigar1);
    }
}
