package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;


import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVDiscoveryTestDataProvider;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;

public class BreakpointComplicationsUnitTest extends GATKBaseTest {

    @Test(groups = "sv")
    public void testGetHomology() {
        final byte[] contigSequence = "ATCGATCGAAAAGCTAGCTA".getBytes();

        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("1", 1, 12), 1, 12, TextCigarCodec.decode("12M8S"), true, 60, 1, 100, ContigAlignmentsModifier.AlnModType.NONE);            // dummy test data, almost guaranteed to be non-factual
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("1", 101, 112), 9, 20, TextCigarCodec.decode("8H12M"), false, 60, 1, 100, ContigAlignmentsModifier.AlnModType.NONE);    // dummy test data, almost guaranteed to be non-factual

        Assert.assertEquals(BreakpointComplications.inferHomology(region1, region2, contigSequence), "AAAA");

        final AlignmentInterval region3 = new AlignmentInterval(new SimpleInterval("1", 1, 8), 1, 8, TextCigarCodec.decode("8M"), true, 60, 1, 100, ContigAlignmentsModifier.AlnModType.NONE);            // dummy test data, almost guaranteed to be non-factual
        final AlignmentInterval region4 = new AlignmentInterval(new SimpleInterval("1", 101, 108), 13, 20, TextCigarCodec.decode("8M"), false, 60, 1, 100, ContigAlignmentsModifier.AlnModType.NONE);    // dummy test data, almost guaranteed to be non-factual

        Assert.assertTrue(BreakpointComplications.inferHomology(region3, region4, contigSequence).isEmpty());
    }

    @Test(groups = "sv")
    public void testGetInsertedSequence() {
        final byte[] contigSequence = "GACGAACGATTTGACTTTAATATGAAATGTTTTATGTGGGCTATAAAATTATCCAAACTCGACACAGGACATTTTGAGCTTATTTCCAAATCATCTGGCCTTCATCTACCCACTGGAACTATTACTCTGCTGGGTCCTCATGGAAACATATCTTTCAGCCCTAACAATGAGACTACAGACATCTACGTCCCCAACACAACAGCTAAAAAGCAGTAGAATGTCAGAAAGGCTATCCACTTAGCCCTTGGCTGACAGGCCCCACTGAGCATCCTTTGCGAAGTCCATTTACTAGCTAATTCATAATTTACACAAGGCATTCAGACATAGCAGCTAAGATATAAAACATTTATCAACACAGGGACTAGTTTGTCATTTTAAAATAATTATGTTTAAGTAAGCCAATAAAGTCTATCTTCTCCAATTTACTTATTGAGCTTTATGAGGCAATTTAAGTCCCGATTTTGGGGGGTATGTATGAAAGGAGAGCATGGAAATGCCATTTGCTCCCTGAAGTTTTTATCTTTTTTTTTTTGAGATAGAGTCTTGTGTTTTCTGTGGAGTACATGAGTATGCATCAAAGCTAACAACGCCCACTGCCCTGTTAGTCAAATACCTTTGA".getBytes();
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("8", 118873207, 118873738), 1, 532, TextCigarCodec.decode("532M87S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("1", 175705642, 175705670), 519, 547, TextCigarCodec.decode("518S29M72S"), false, 3, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region3 = new AlignmentInterval(new SimpleInterval("1", 118875262, 118875337), 544, 619, TextCigarCodec.decode("543S76M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);

        Assert.assertTrue(BreakpointComplications.inferInsertedSequence(region3, region1, contigSequence).isEmpty());
        Assert.assertEquals(BreakpointComplications.inferInsertedSequence(region1, region3, contigSequence), "GAGATAGAGTC");

        Assert.assertTrue(BreakpointComplications.inferInsertedSequence(region2, region1, contigSequence).isEmpty() && BreakpointComplications.inferInsertedSequence(region1, region2, contigSequence).isEmpty());
    }

    // -----------------------------------------------------------------------------------------------
    // Tests for CIGAR extraction on tandem duplications
    // -----------------------------------------------------------------------------------------------
    @Test(groups = "sv")
    public void testExtractCigarForSimpleTandup() {

        final int contigTotalLength = 355;

        // forward strand
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("1", 1000001, 1000125), 16, 90,
                TextCigarCodec.decode("5H10S15M20D25M30D35M260S5H"),
                true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("1", 1000041, 1000145), 191, 345,
                TextCigarCodec.decode("5H185S45M30I55M20I5M10S5H"),
                true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);


        final Cigar cigar1 = BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications.extractCigarForTandupExpansion(region1, 1000125, 1000041);
        Assert.assertEquals(cigar1, TextCigarCodec.decode("20M30D35M"));
        final Cigar cigar2 = BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications.extractCigarForTandupExpansion(region2, 1000125, 1000041);
        Assert.assertEquals(cigar2, TextCigarCodec.decode("45M30I40M"));

        // reverse strand
        final AlignmentInterval region3 = new AlignmentInterval(region2.referenceSpan, contigTotalLength-region2.endInAssembledContig+1, contigTotalLength-region2.startInAssembledContig+1,
                CigarUtils.invertCigar(region2.cigarAlong5to3DirectionOfContig),
                false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region4 = new AlignmentInterval(region1.referenceSpan, contigTotalLength-region1.endInAssembledContig+1, contigTotalLength-region1.startInAssembledContig+1,
                CigarUtils.invertCigar(region1.cigarAlong5to3DirectionOfContig),
                false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);

        final Cigar cigar3 = BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications.extractCigarForTandupExpansion(region3, 1000125, 1000041);
        Assert.assertEquals(CigarUtils.invertCigar(cigar3), cigar2);
        final Cigar cigar4 = BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications.extractCigarForTandupExpansion(region4, 1000125, 1000041);
        Assert.assertEquals(CigarUtils.invertCigar(cigar4), cigar1);

        // modifying the real event below
        //  "asm030282:tig00005     1_342_chrX:1294785-1295169_-_168M43D174M1321H_60_65_173_O       463_1663_chrX:1293941-1295139_-_462S1098M2I101M_60_17_1106_O";
        final AlignmentInterval region5 = new AlignmentInterval(new SimpleInterval("chr20", 1294785, 1295169),
                1, 342, TextCigarCodec.decode("168M43D174M1321H"),
                false, 60, 65, 173, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region6 = new AlignmentInterval(new SimpleInterval("chr20", 1293941, 1295139),
                463, 1663, TextCigarCodec.decode("462S1098M2I101M"),
                false, 60, 17, 1106, ContigAlignmentsModifier.AlnModType.NONE);
        final SimpleChimera simpleChimera = new SimpleChimera(region5, region6, Collections.emptyList(), "asm030282:tig00005",
                SimpleSVDiscoveryTestDataProvider.b38_seqDict);
        final BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications smallDuplicationWithPreciseDupRangeBreakpointComplications = new BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications(simpleChimera, "CTCCTGTCGTCCAGGTAGACATGGAGCAGGCTCTCCTTCTTGTCTACACTGGGTCCAGGTGGTTGTGGGACAAGATCTCCTCTTGTCTACACTGGATCGAGGTGGATGTGAGGCAGTCTCTCTTCCTGTCTACACTGGGTCCAGGTGGTTGTGGGACAAGAGCTCCTCCTGTCTACACTGGGTCCAGGTAGCTGTGGGACAAGAGCTCCTCCTGTCTACACTGGGTCTCCATGGAGGTGGGGCAGGGTCTCCTTCTGTCTACACTGCGTGTAGTTGGAGGTGGGGCAGGGTCTCCTCCTGTCTACACTGGGTCCAGGTAGACATGGGGCAGTCTCTCCTTCTTGTCTACACTGGGTCCAGGTGGTTGTGGGACAAGAGCTCCTCCTGTCTACACTGGGTCCAGGTAGTTGTGGGACAAGAGCTCCTCCTGTCTACACTGGGTCTCCATGGAGGTGGGGCAGTCTCTCCTTCTTGTCTACACTGGGTCCAGGTGGTTGTGGGACAAGATCTCCTCTTGTCTACACTGGCTCGAGGTGGACATGGGGCAGGGTCTCTTCTTGTCTACACTGGGTCCAGGAGGTTGTGAGACAAGATCTCCTCTTGTCTACACTGGATCGAGGTGGACGTGAGGCAGTCTCTCTCCCTGTCTACACTGGGTCCAGGTAGTTGTGGGACAAGAGCTCCTCCTGTCTACACTGGGTCTCCATGGAGGTGGAGCAGGGTCTCCTCCTGTCTACACTGGGTGTAGGTGGAGGTGGGGTCGGGTGTCCTCCTATCTACACTGGGTCCAGGTAGACATGGGGCAGGGTCTCCTTCTCTCTACACTGCGTCCAGCTGGAGGTGGAGCAGAGGCTCTCCTTGCTTGTGGCATCGTCCCCCACACCTCCCGGTCCACTTCCTGGTTCCATGGTTGCAGGATCATCCTTGTCCACCCTCCCTGCAACCTCTTTCAAGGTGGCTCCACAGGCCACAGACCCTTCACCTCTTCCTCCGCTACCCGAAGTGTGTTCACCCCAGAGTCACCGCTCACCACCCACCCATCCTTCCCCCAGGCCACTTCCCCGGGATTCCCAGGCTCCTGTGCGGGCGTGTCCCGTACGCCTCCCTCCTGGTGCCCAGCCCCGGGGAGCTCTCACCGACCTTTCTGTGGACGTCCAGCTGGTACTGAAAGTCCAGGTACGACAGCTTCTGATAGGTCCTGGGCTGTTTCCACCGTACGAGGCAGTGCGTCGTGTTGCAACGTACGGTGACATTGCTGGGAGGGTTGAATCGTTCTGTAACGAGGGCGCAGGACACACCCCTGAACCCGAGAGGTCCTGTCTACACTGGGTCCAGGTGGAGGTGGTGCAGAGTCTCCTCCTGTCTACACTGGGTCCAGGTGGAGGTGGAGTAGGTCCTGTCTACACTGGGTCCAGGTGGAGGTGGAGTAGGGACACCTCTTTGACTACACTGGGTCCAGGTGGAGATGGGGCAGGGTCTCCTCCTGTCTACACTGGGTCCAGGTGGAGGTGGGGCAGGGTCTCCTCCTGTCTACACTGGGTCTAGGTGGAGGTGGTGCAGAGTCTTCTCCTGTCTACACTGGGTCCAGGTGGAGGTGGGGCAGGGTCTCCTCCTGTCTACACTCGGTCCAGGTGGATGTGGACTAGGGACACCTCTTTGTCTA".getBytes());
        Assert.assertEquals(smallDuplicationWithPreciseDupRangeBreakpointComplications.getCigarStringsForDupSeqOnCtgForwardStrandRep(),
                Arrays.asList("355M", "174M43D138M"));
    }
}
