package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.*;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;


public class AssemblyAlignmentParserUnitTest extends BaseTest{
    private static final PipelineOptions dummyOptions = null;
    private static final SAMSequenceDictionary seqDict = new ReferenceMultiSource(dummyOptions, b37_reference_20_21, ReferenceWindowFunctions.IDENTITY_FUNCTION).getReferenceSequenceDictionary(null);

    private static final String chrForLongContig1 = "20";


    @Test
    public void testGetInsertedSequence() {
        final byte[] contigSequence = "GACGAACGATTTGACTTTAATATGAAATGTTTTATGTGGGCTATAAAATTATCCAAACTCGACACAGGACATTTTGAGCTTATTTCCAAATCATCTGGCCTTCATCTACCCACTGGAACTATTACTCTGCTGGGTCCTCATGGAAACATATCTTTCAGCCCTAACAATGAGACTACAGACATCTACGTCCCCAACACAACAGCTAAAAAGCAGTAGAATGTCAGAAAGGCTATCCACTTAGCCCTTGGCTGACAGGCCCCACTGAGCATCCTTTGCGAAGTCCATTTACTAGCTAATTCATAATTTACACAAGGCATTCAGACATAGCAGCTAAGATATAAAACATTTATCAACACAGGGACTAGTTTGTCATTTTAAAATAATTATGTTTAAGTAAGCCAATAAAGTCTATCTTCTCCAATTTACTTATTGAGCTTTATGAGGCAATTTAAGTCCCGATTTTGGGGGGTATGTATGAAAGGAGAGCATGGAAATGCCATTTGCTCCCTGAAGTTTTTATCTTTTTTTTTTTGAGATAGAGTCTTGTGTTTTCTGTGGAGTACATGAGTATGCATCAAAGCTAACAACGCCCACTGCCCTGTTAGTCAAATACCTTTGA".getBytes();
        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("8", 118873207, 118873739), TextCigarCodec.decode("532M87S"), true, 60, 0, 1, 532);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 175705642, 175705671), TextCigarCodec.decode("518S29M72S"), false, 3, 0, 519, 547);
        final AlignmentRegion region3 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 118875262, 118875338), TextCigarCodec.decode("543S76M"), false, 60, 0, 544, 619);

        Assert.assertTrue(AssemblyAlignmentParser.getInsertedSequence(region3, region1, contigSequence).isEmpty());
        Assert.assertEquals(AssemblyAlignmentParser.getInsertedSequence(region1, region3, contigSequence), "GAGATAGAGTC");

        Assert.assertTrue(AssemblyAlignmentParser.getInsertedSequence(region2, region1, contigSequence).isEmpty() && AssemblyAlignmentParser.getInsertedSequence(region1, region2, contigSequence).isEmpty());
    }

    @Test
    public void testARtooSmall() {
        final byte[] contigSequence = SVVariantConsensusCallUnitTest.LONG_CONTIG1.getBytes();
        AlignmentRegion region1 = new AlignmentRegion("702700", "702700", new SimpleInterval(chrForLongContig1, 20138007, 20142231), TextCigarCodec.decode("1986S236M2D1572M1I798M5D730M1I347M4I535M"), false, 60, 36, 1, contigSequence.length - 1986);
        AlignmentRegion region2 = new AlignmentRegion("702700", "702700", new SimpleInterval(chrForLongContig1, 20152030, 20154634), TextCigarCodec.decode("3603H24M1I611M1I1970M"), true, 60, 36, 3604, contigSequence.length);
        Assert.assertFalse( AssemblyAlignmentParser.currentAlignmentIsTooShort(region1, region2, SVConstants.DEFAULT_MIN_ALIGNMENT_LENGTH) );
        Assert.assertFalse( AssemblyAlignmentParser.currentAlignmentIsTooShort(region2, region1, SVConstants.DEFAULT_MIN_ALIGNMENT_LENGTH) );

        Assert.assertFalse( AssemblyAlignmentParser.currentAlignmentIsTooShort(region1, region2, 3000) );
        Assert.assertTrue( AssemblyAlignmentParser.currentAlignmentIsTooShort(region2, region1, 3000) );
    }

    @Test
    public void testGetAssembledBreakpointFromAlignmentRegionsStrangeLeftBreakpoint() throws Exception {
        final byte[] contigSequence = SVVariantConsensusCallUnitTest.LONG_CONTIG1.getBytes();
        AlignmentRegion region1 = new AlignmentRegion("702700", "702700", new SimpleInterval(chrForLongContig1, 20138007, 20142231), TextCigarCodec.decode("1986S236M2D1572M1I798M5D730M1I347M4I535M"), false, 60, 36, 1, contigSequence.length - 1986);
        AlignmentRegion region2 = new AlignmentRegion("702700", "702700", new SimpleInterval(chrForLongContig1, 20152030, 20154634), TextCigarCodec.decode("3603H24M1I611M1I1970M"), true, 60, 36, 3604, contigSequence.length);
        final List<AlignmentRegion> alignmentRegionList = Arrays.asList(region1, region2);
        final List<ChimericAlignment> assembledBreakpointsFromAlignmentRegions = AssemblyAlignmentParser.getChimericAlignmentsFromAlignmentRegions(new Tuple2<>(alignmentRegionList, contigSequence));
        Assert.assertEquals(assembledBreakpointsFromAlignmentRegions.size(), 1);
        final ChimericAlignment chimericAlignment = assembledBreakpointsFromAlignmentRegions.get(0);
        Assert.assertEquals(chimericAlignment.contigId, "702700");
        Assert.assertEquals(chimericAlignment.regionWithLowerCoordOnContig, region1);
        Assert.assertEquals(chimericAlignment.regionWithHigherCoordOnContig, region2);
        Assert.assertFalse(chimericAlignment.homology.isEmpty());

        final Tuple2<SimpleInterval, SimpleInterval> leftAndRightBreakpointsOnReferenceLeftAlignedForHomology = chimericAlignment.getLeftJustifiedBreakpoints(seqDict);

        Assert.assertEquals(leftAndRightBreakpointsOnReferenceLeftAlignedForHomology._1(), new SimpleInterval(chrForLongContig1, 20138007, 20138007));
        Assert.assertEquals(leftAndRightBreakpointsOnReferenceLeftAlignedForHomology._2(), new SimpleInterval(chrForLongContig1, 20152651, 20152651));
    }

    @Test
    public void testGetAssembledBreakpointsFromAlignmentRegions() throws Exception {
        final byte[] contigSequence = "GACGAACGATTTGACTTTAATATGAAATGTTTTATGTGGGCTATAAAATTATCCAAACTCGACACAGGACATTTTGAGCTTATTTCCAAATCATCTGGCCTTCATCTACCCACTGGAACTATTACTCTGCTGGGTCCTCATGGAAACATATCTTTCAGCCCTAACAATGAGACTACAGACATCTACGTCCCCAACACAACAGCTAAAAAGCAGTAGAATGTCAGAAAGGCTATCCACTTAGCCCTTGGCTGACAGGCCCCACTGAGCATCCTTTGCGAAGTCCATTTACTAGCTAATTCATAATTTACACAAGGCATTCAGACATAGCAGCTAAGATATAAAACATTTATCAACACAGGGACTAGTTTGTCATTTTAAAATAATTATGTTTAAGTAAGCCAATAAAGTCTATCTTCTCCAATTTACTTATTGAGCTTTATGAGGCAATTTAAGTCCCGATTTTGGGGGGTATGTATGAAAGGAGAGCATGGAAATGCCATTTGCTCCCTGAAGTTTTTATCTTTTTTTTTTTGAGATAGAGTCTTGTGTTTTCTGTGGAGTACATGAGTATGCATCAAAGCTAACAACGCCCACTGCCCTGTTAGTCAAATACCTTTGA".getBytes();
        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("8", 118873207, 118873739), TextCigarCodec.decode("532M87S"), true, 60, 0, 1, 532);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 175705642, 175705671), TextCigarCodec.decode("518S29M72S"), false, 3, 0, 519, 547);
        final AlignmentRegion region3 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 118875262, 118875338), TextCigarCodec.decode("543S76M"), false, 60, 0, 544, 619);
        final List<AlignmentRegion> alignmentRegionList = Arrays.asList(region1, region2, region3);
        final List<ChimericAlignment> assembledBreakpointsFromAlignmentRegions = AssemblyAlignmentParser.getChimericAlignmentsFromAlignmentRegions(new Tuple2<>(alignmentRegionList, contigSequence));
        Assert.assertEquals(assembledBreakpointsFromAlignmentRegions.size(), 1);
        final ChimericAlignment chimericAlignment = assembledBreakpointsFromAlignmentRegions.get(0);
        Assert.assertEquals(chimericAlignment.contigId, "contig-1");
        Assert.assertEquals(chimericAlignment.regionWithLowerCoordOnContig, region1);
        Assert.assertEquals(chimericAlignment.regionWithHigherCoordOnContig, region3);
        Assert.assertEquals(chimericAlignment.homology, "");
        Assert.assertEquals(chimericAlignment.insertedSequence, "GAGATAGAGTC");
    }

    @Test
    public void testGetAssembledBreakpointsFromAlignmentRegionsWithOverlappingAlignmentRegion() throws Exception {
        final byte[] contigSequence = "ACTAGAGCATCTACGTGTTCCTGTGGTTTTGGAGCAAGAGTGATTTGAGTTTCAGAGATTTTTACTAATTCTTCTTCCCCTACCAGAAAAAAAGATCTTACCATTTGAGAGTGAGATGTAAACCCAGCCCTGTCTGACCTGAGTCTGTGCCCTAAGCCTATGCTAAGCCAAGCAGTGCCTGGAGCCACCACAGGTCCACACAATTCGTTAACATGATGAAGCAAGGATGGAAATTGGACAAAATAGTGTGCCTACTGAATCTAAGAATGAAAAATGATTGCACTCCTACTCTGAGTGCTTTGGAGCACTGCCCAGTTGGGCAAAGGGTCAGCGCCTGGGCAGAGGTCCCCACAACCTGGCAGGAGTGTGGTCGGCCACCCTATGGGCCTCCATCATGTGCAGTGACAGCGGGGCTGTCATGTCACCGTGTGGGAGGGCTTGCAGGTGAAGTGGTCTGGGAGGGGTCCCCCAGACAAAGCCAAGGTTCTGAGAGTTGGCCCGAACACTGCTGGATTCCACTTCACCTGCAAGCCCTCCCACACGGTGACATGACAGCCTATAATACAGTTCCGCATGGCCACGTCATACAACCCTGTCATATTGGTGAGCAATTGCTGTGTAGCCAAAGACCCCAAAACTCAAACAGCATTTATTATTATTGCCCCCATGTCTGAGAGTCAGATGTGCATTTGCTGATCTCAGCTTGTTTGAGCTGCTGCAGGGTTGGGGCTCTGCTCCAGGCAGGCTTAGCTGTCACCACATGCACACATACATTCTGGGCCTCTGCTGCGCGCGTCACGTTCACTGAAGATCTTGGGATTGGGAGTTAGGGCGGTGGGAGGGCCCAGCAAAGTCACCTGGCGATGGCAGGGACACAGGGAGGAATGTAGAATGGGGCCGATGATGGGACCCACACGTCTGCAAAGCTGCGGTCTCCTTGAGGGGTGGAGACAGCAACAACTCACCGCACGCGGTGCTTCAGTTCACCATCTCCCTGGGACATTAGGGGGCCCCGTGTTATCTCATTTTGCTCTGGTTTGCATTAGTTTTTTATCACTTCGTAGATGAAGCCACTGACACCCAGAGAGGGAAAGTGGCCTGACCAAGGGCCACAGCAGGGGAGCGAAGGAGCCCCACAGTTCGGCAGGAACACAGCCTCTCCCTGGCTTTCAGGTTCACTGACATCTTCTCATGGCCTCTGTAACTCACCAGGCATCAGGGTGTAGTCCTTAGACCAGTGTCCCACAGCTGCCACAGAGTGGGAGCTCACCATCAGTTATAAGTCACTAGAAAGGCTTTTGGACATTATAAGCTACAATGGAAAATAAGTCATCTGTGGATTTTTGTGACAGATTCCAAAAATTTGAATATTTTGTCTACTTAGGTTTTTGGTTAATTTTATCCTCAAAACTGTTCTGCAGTGATTAAGCTGTACAAACTGCATCATGGGCGAATTGGCATATTCAGAAATGACTGATATTCTTGATTTCAGTTTTTTACTTTGTATGTAGCTCCTCAAGGAAAC".getBytes();
        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 23102817, 23103304), TextCigarCodec.decode("487M1006S"), true, 60, 1, 1, 487);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 23103196, 23103238), TextCigarCodec.decode("483S42M968S"), false, 60, 2, 484, 525);
        final AlignmentRegion region3 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 23103633, 23104603), TextCigarCodec.decode("523S970M"), true, 60, 3, 524, 1493);
        final List<AlignmentRegion> alignmentRegionList = Arrays.asList(region1, region2, region3);
        final List<ChimericAlignment> assembledBreakpointsFromAlignmentRegions = AssemblyAlignmentParser.getChimericAlignmentsFromAlignmentRegions(new Tuple2<>(alignmentRegionList, contigSequence));
        Assert.assertEquals(assembledBreakpointsFromAlignmentRegions.size(), 1);
        final ChimericAlignment chimericAlignment = assembledBreakpointsFromAlignmentRegions.get(0);
        Assert.assertEquals(chimericAlignment.contigId, "contig-1");
        Assert.assertEquals(chimericAlignment.regionWithLowerCoordOnContig, region1);
        Assert.assertEquals(chimericAlignment.regionWithHigherCoordOnContig, region3);
        Assert.assertEquals(chimericAlignment.homology, "");
        Assert.assertEquals(chimericAlignment.insertedSequence, "TGAGAGTTGGCCCGAACACTGCTGGATTCCACTTCA");
        Assert.assertEquals(chimericAlignment.insertionMappings.size(), 1);
        Assert.assertEquals(chimericAlignment.insertionMappings.get(0), "1-contig-1:484-525:20,23103196,-,483S42M968S,60,2");
    }

    @Test
    public void testTreatAlignmentRegionAsInsertion() throws Exception {
        AlignmentRegion overlappingRegion1 = new AlignmentRegion("overlap", "22", new SimpleInterval("19", 48699881, 48700035), TextCigarCodec.decode("47S154M"), false, 60, 0, 1, 154);
        AlignmentRegion overlappingRegion2 = new AlignmentRegion("overlap", "22", new SimpleInterval("19", 48700584, 48700669), TextCigarCodec.decode("116H85M"), true, 60, 0, 117, 201);

        Assert.assertTrue(AssemblyAlignmentParser.nextAlignmentMayBeNovelInsertion(overlappingRegion1, overlappingRegion2, 50));
    }

    @Test
    public void testGetHomology() {
        final byte[] contigSequence = "ATCGATCGAAAAGCTAGCTA".getBytes();

        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 1, 12), TextCigarCodec.decode("12M8S"), true, 60, 1, 1, 12);            // dummy test data, almost guaranteed to be non-factual
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 101, 112), TextCigarCodec.decode("8H12M"), false, 60, 1, 9, 20);    // dummy test data, almost guaranteed to be non-factual

        Assert.assertEquals(AssemblyAlignmentParser.getHomology(region1, region2, contigSequence), "AAAA");

        final AlignmentRegion region3 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 1, 12), TextCigarCodec.decode("8M"), true, 60, 1, 1, 8);            // dummy test data, almost guaranteed to be non-factual
        final AlignmentRegion region4 = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 101, 112), TextCigarCodec.decode("8M"), false, 60, 1, 13, 20);    // dummy test data, almost guaranteed to be non-factual

        Assert.assertTrue(AssemblyAlignmentParser.getHomology(region3, region4, contigSequence).isEmpty());
    }

    @Test
    public void testGappedAlignmentBreaker_OneInsertion() {

        final Cigar cigar = TextCigarCodec.decode("56S27M15I32M21S");
        final AlignmentRegion alignmentRegion = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 158), cigar, true, 60, 0, 57, 130);

        final List<AlignmentRegion> generatedARList = StreamSupport.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, 1).spliterator(), false).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.size(), 2);
        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 126), TextCigarCodec.decode("56S27M68S"), true, -60, -1, 57, 83));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 127, 158), TextCigarCodec.decode("98S32M21S"), true, -60, -1, 99, 130));
    }

    @Test
    public void testGappedAlignmentBreaker_OneDeletion() {
        final Cigar cigar = TextCigarCodec.decode("2S205M2D269M77S");
        final AlignmentRegion alignmentRegion = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 575), cigar, true, 60, 0, 208, 476);

        final List<AlignmentRegion> generatedARList = StreamSupport.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, 1).spliterator(), false).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.size(), 2);
        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 304), TextCigarCodec.decode("2S205M346S"), true, -60, -1, 3, 207));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 307, 575), TextCigarCodec.decode("207S269M77S"), true, -60, -1, 208, 476));
    }

    @Test
    public void testGappedAlignmentBreaker_Complex() {

        final AlignmentRegion alignmentRegion = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 414), TextCigarCodec.decode("397S118M2D26M6I50M7I26M1I8M13D72M398S"), true, 60, 65, 398, 711);

        final List<AlignmentRegion> generatedARList = StreamSupport.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, 1).spliterator(), false).collect(Collectors.toList());

        Assert.assertEquals(generatedARList.size(), 6);

        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 217), TextCigarCodec.decode("397S118M594S"), true, -60, -1, 398, 515));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 220, 245), TextCigarCodec.decode("515S26M568S"), true, -60, -1, 516, 541));
        Assert.assertEquals(generatedARList.get(2), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 246, 295), TextCigarCodec.decode("547S50M512S"), true, -60, -1, 548, 597));
        Assert.assertEquals(generatedARList.get(3), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 296, 321), TextCigarCodec.decode("604S26M479S"), true, -60, -1, 605, 630));
        Assert.assertEquals(generatedARList.get(4), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 322, 329), TextCigarCodec.decode("631S8M470S"), true, -60, -1, 632, 639));
        Assert.assertEquals(generatedARList.get(5), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 343, 414), TextCigarCodec.decode("639S72M398S"), true, -60, -1, 640, 711));
    }

    @Test
    public void testGappedAlignmentBreaker_Sensitivity() {

        final Cigar cigar = TextCigarCodec.decode("10M10D10M60I10M10I10M50D10M");
        final AlignmentRegion alignmentRegion = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 209), cigar, true, 60, 0, 1, 120);

        final List<AlignmentRegion> generatedARList = StreamSupport.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, SVConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY).spliterator(), false).collect(Collectors.toList());

        Assert.assertEquals(generatedARList.size(), 3);
        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 129), TextCigarCodec.decode("10M10D10M100S"), true, -60, -1, 1, 20));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 130, 149), TextCigarCodec.decode("80S10M10I10M10S"), true, -60, -1, 81, 110));
        Assert.assertEquals(generatedARList.get(2), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 200, 209), TextCigarCodec.decode("110S10M"), true, -60, -1, 111, 120));
    }

    @Test
    public void testGappedAlignmentBreaker_HardAndSoftClip() {

        final AlignmentRegion alignmentRegion = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 138), TextCigarCodec.decode("1H2S3M5I10M20D6M7S8H"), true, 60, 0, 4, 27);

        final List<AlignmentRegion> generatedARList = StreamSupport.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, 1).spliterator(), false).collect(Collectors.toList());

        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 102), TextCigarCodec.decode("1H2S3M28S8H"), true, -60, -1, 4, 6));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 103, 112), TextCigarCodec.decode("1H10S10M13S8H"), true, -60, -1, 12, 21));
        Assert.assertEquals(generatedARList.get(2), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 133, 138), TextCigarCodec.decode("1H20S6M7S8H"), true, -60, -1, 22, 27));
    }
}
