package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.ContigAligner.AlignmentRegion;
import static org.broadinstitute.hellbender.tools.spark.sv.ContigAligner.AssembledBreakpoint;

public class ContigAlignerTest extends BaseTest {

    private ContigAligner contigAligner;

    @BeforeClass
    public void setup() throws Exception {
        contigAligner = new ContigAligner(b37_reference_20_21);
    }

    @Test
    public void testAlignContigs() throws Exception {

        // data taken from G94982 NA12878 assembly for breakpoint interval 21:27373582-27375176, which contains a small inversion
        final List<String> contigsData = new ArrayList<>();
        contigsData.add(">contig-16 169 0");
        contigsData.add("ACATGTTAATTAAGTAAATAGTATAAGCTACACGTTGAGCCTACAACTCTTGCTGGAGTCTCCAGCTGTAAGAGCCTTAACTCCTTCCTTTCATAAAGGCCAAGAAAAGAGTGAGTTAGTCCCTATTCGGGGAGTGTTGGTTGGTTCTCCCAATCCTGTTAGGGGCAGT");
        contigsData.add(">contig-5 151 0");
        contigsData.add("GAGTCTATATATTAAATAAAAGTGGTAGAATTTTGGATGGGTGCAAGTCAAGTAATGGGATAAAAGAATGATGATCGAAGGAGATGATTCATGATAACTATGTGCGCTTATGTGAAACCTAGCAACAATTTCTTACTTATGATAAAGGCAC");
        contigsData.add(">contig-11 192 0");
        contigsData.add("GTGAACATCACAGCTCTGCAGATGCACAAACCTTGTCTTTATGAGTCTATATATTAAATAAAAGTGGTAGAATTTTGGATGGGTGCAAGTCAAGTAATGGGATAAAAGAATGATGATCGAAGGAGATGATTCATGATAACTATGTGCCCTTATGTGAAACCTAGCAACAATTTCTTACTTATGATAAAGGCA");
        contigsData.add(">contig-9 2009 0");
        contigsData.add("CATGAAGTACATCTATAACAACGTATTATTAATTCTGCTATTTAATTGATATTTATATATTAAAAATCTGTATTATTTGTATATATTTACATACACATTAAGTACAACTATAACAATGTATTATTTAATTATGCTACTTAATTGATATTTATGTATTAAAAATCTGTATTATTTGTATATTTTGTATATTTACCTGACCCTTGGGCAGGAATGTTCCAAACAATAAAACGCAGTCAATCAACTCCAACTGATTAATACACACCACGCACTGTTATACATTAAAACACATCTTACCTGTCAACATGTCTTGCCACCCTCCCCAAGATCAGTTTCACTCATCCCTGATTGCCTCAAAGTCAATCAGTGACCCATCCAATTGCATGTATAAAGTCAATCCTGAGTAGAATTATAACTCTACCTCCCATGTATATCCCATTTTTATTGGAAGTGAAACAACAGCCCAATTCTCTCTCTGGGTTTGTGCCCTGCCCTTGCTAGACTGACGCCTTAGCACACACAGATGTCTTCTTAGTTTTACTGACAGAGGAAAAAACACCTTTCTAAAAGATGAAAGCCACGCATCACCCATTACTTCCAAGTTAGGAATATTCAACAGATTCTTTTTTTCTTTTTATGAGACAGGGTCTTGCACTGCCGCCCAGTGGCACAATCTTGGCTCACTGCAACCTCTACCTCCTGAGTTCAAGTGATTCTTGCGCCTCAGCCTCCCAAGTAGCTGGGATTACCCTTGTACCAGCACGCCTGGCTAATTTTTTGGTATTTTTAGTAGAGATGGGGTTTCACCAGGTAGGCCAGGCTGGTCTTGAACTCCTGGCCTCAAGTGATACACCTGCCTCAGCCTCTCTAAGTGCTGGGATTACAGGCGTGAGCCACTGCCACTGGTCCTTCAACAGATTGATTCTAATTAGCCAATCAAAGACAAGGATCCACAATGTCAATACAAATGGGGACTAACTATTGTTTTACTTCCCTATACACGTTGCTCTTTCAGTAGATGCAATAAAGTACTGGTAAAACCAGAGGTGGCTACCATCACGATGATGTCAACAGGAGGGACAGTCAGCACTAAGCCCAGAAGGTGTCAAACACTCCGCAGGAGAAATGCGCCATGCAACGGACATGAAGATGATCTGACACTCTTCACGTGGTTTTCAGATGGAAACGTGGCTACGAAAGCATCAACCTCATTATCATCCATCATTAAGGCCATCTCACTCAGTACTGCTGCTTTCAAAGTCCACGCTCCCAAAGCAAATTGGATTTCTGTACACAATACTCTTACAGGATGAAACCCAACCAACTTGTTGACGTAAGTATCCCATCTACTTCGCAATTTTATTTATCTTCCAAAATTAAAGGACTGGCACCCTGATTTATTAAAAGTGAATTGGTTCTAGGGACCATATCCCCTCTGAGTTACTGACAGAGCAGCTTCTGGCCTGTGAAGCTCAAAGCCATGCCTAGATGTGAGGATCCCATCAACTATAGCCAGCAGTGGTCTCTGCTCCCAGAGCTTACCTTCAGTTGTAACAAACATCTACAGCTTCACTTTGCTTTCCTTCTCTGTATCTCCTTCACCAATGTGTACATATTCATGTCACATGCTCTTTAACGTTTCTAAGACACAACAATACCCATTATTAATCTCTTACTCAGGCAGTGTTAATTCCAATAAGACTGCAAGGCAGCAGTGCCTCCAGAGCCTGCACTGTGCCGGGAATTGGCACTGGGATTGCAAGCACCGTGGTCATTGCTACCAAGGGAGGCACAGAATCCCTTCACCCCATAGTCACGGGAGCATTGGCAACAAAATTTACCATGCCTTGGCAATTTATGCTGACCCTCACATTTTGGCATTAAAAAAAAGTATCATCTTAGAGTACCAAATATGTGTCTTCTTATAGATGTGTGCCACGGGCCTCGCAGTGAACATCACAGCTCTGCAGATGCACAAACCTTGTCTTTATGAGTCTATATATTAAATAAAA");
        contigsData.add(">contig-2 152 0");
        contigsData.add("TCACATTTTGGCATTAAAAAAAAGTATCATCTTAGAGTACCAAATATGTGTCTTCTTATAGATGTGTGCCACGGGCCTCGCAGTGAACATCACAGCTCTGCAGATGCACAAACCTTGTCTTTATGAGTCTATATATTAAATAAAAGTGGTAG");
        contigsData.add(">contig-13 165 0");
        contigsData.add("GTGTGCCACGGGCCTCGCAGTGAACATCACAGCTCTGCAGATGCACAAACCTTGTCTTTATGAGTCTATATATTAAATAAAAGTGGTAGAATTTTGGATGGGTGCAAGTCAAGTAATGGGATAAAAGAATGATGATCGAAGGAGATGATTCATGATAACTATGTG");
        contigsData.add(">contig-0 250 0");
        contigsData.add("CATTCTATAACAGCTATGAAACCATTTGTTATAGCTAATTCGTAGAAGCCAGGCTAAATGAGCCCTAATCAACTAGAACAAGCATTCAGAACTCAAAATCCATTTACATGTTAATTAAGTAAATAGTATAAGCTACACGTTGAGCCTACAACTCTTGCTGGAGTCTCCAGCTGTAAGAGCCTTAACTCCTTCCTTTCATAAAGGCCAAGAAAAGAGTGAGTTAGTCCCTATTCGGGGAGTGTTGGTTGGT");
        contigsData.add(">contig-15 151 0");
        contigsData.add("ACAAACCTTGTCTTTATGAGTCTATATATTAAATAAAAGTGGTAGAATTTTGGATGCGTGCAAGTCAAGTAATGGGATAAAAGAATGATGATCGAAGGAGATGATTCATGATAACTATGTGCCCTTATGTGAAACCTAGCAACAATTTCTT");
        contigsData.add(">contig-6 241 0");
        contigsData.add("CTTATGTGAAACCTAGCAACAATTTCTTACTTATGATAAAGGCACAATCTCCCCAGTCCAGAATCTAAATGGCATTATATTAACATTTTCCTATACCACTAATTAAACATAATGCTACATTCATTTTATTTTCTTCTGCCTCATTCTATAACAGCTATGAAACCATTTGTTATAGCTAATTCGTAGAAGCCAGGCTAAATGAGCCCTAATCAACTAGAACAAGCATTCAGAACTCAAAATC");
        contigsData.add(">contig-10 151 0");
        contigsData.add("TGTGTCTTCTTATAGATGTGTGCCACGGGCCTCGCAGTGAACATCACAGCTCTGCAGATGCACAAACCTTGTCTTTATGAGTCTATATATTAAATAAAAGAGGTAGAATTTTGGATGGGTGCAAGTCAAGTAATGGGATAAAAGAATGATG");
        contigsData.add(">contig-17 151 0");
        contigsData.add("TAGAAGCCAGGCTAAATGAGCCCTAATCAACTAGAACAAGCATTCAGAACTCAAAATCCCTTTACATGTTAATTAAGTAAATAGTATAAGCTACACGTTGAGCCTACAACTCTTGCTGGAGTCTCCAGCTGTAAGAGCCTTAACTCCTTCC");
        contigsData.add(">contig-3 195 0");
        contigsData.add("AAGTCAAGTAATGGGATAAAAGAATGATGATCGAAGGAGATGATTCATGATAACTATGTGCCCTTATGTGAAACCTAGCAACAATTTCTTACTTATGATAAAGGCACAATCTCCCCAGTCCAGAATCTAAATGGCATTATATTAACATTTTCCTATACCACTAATTAAACATAATGCTACATTCATTTTATTTTC");
        contigsData.add(">contig-4 151 0");
        contigsData.add("CCTTATGTGAAACCTAGCAACAATTTCTTACTTATGATAAAGGCACAATCTCCCCAGTCCAGAATCTAAATGGCATTATATTAACATTTTCCTATACCACTAATTAAACATAATGCTACATTCATTTTATTTTCTTCTGCCTCATTCTATA");
        contigsData.add(">contig-12 195 0");
        contigsData.add("AAAAAAAGTATCATCTTAGAGTACCAAATATGTGTCTTCTTATAGATGTGTGCCACGGGCCTCGCAGTGAACATCACAGCTCTGCAGATGCACAAACCTTGTCTTTATGAGTCTATATATTAAATAAAAGTGGTAGAATTTTGGATGGGTGCAAGTCAAGTAATGGGATAAAAGAATGATGATCGAAGGAGATGA");
        contigsData.add(">contig-8 194 0");
        contigsData.add("CCACTAATTAAACATAATGCTACATTCATTTTATTTTCTTCTGCCTCATTCTATAACAGCTATGAAACCATTTGTTATAGCTAATTCGTAGAAGCCAGGCTAAATGAGCCCTAATCAACTAGAACAAGCATTCAGAACTCAAAATCCATTTACATGTTAATTAAGTAAATAGTATAAGCTACACGTTGAGCCTA");
        contigsData.add(">contig-1 151 0");
        contigsData.add("AGGCACAATCTCCCCAGTCCAGAATCTAAATGGCATTATATTAACATTTTCCTATACCACTAATTAAACATAATGCTACATTCATTTTATTTTCTTCTGTCTCATTCTATAACAGCTATGAAACCATTTGTTATAGCTAATTCGTAGAAGC");
        contigsData.add(">contig-14 151 0");
        contigsData.add("TAGAATTTTGGATGGGTGCAAGTCAAGTAATGGGATAAAAGAATGATGATCGAAGGAGATGATTCATGATAACTATGTGCCCTTATGTGAAACCTAGCAACAATTTCTTACTTATGATAAAGGCACAATCTCCCCAGTCCAGAATCTAAAT");
        contigsData.add(">contig-7 151 0");
        contigsData.add("TATTAAATAAAAGTGGTAGAATTTTGGATGGGTGCAAGTCAAGTAATGGGATAAAAGAATGATGATCGAAGGAGATGACTCATGATAACTAAGTGCCCTTATGTGAAACCTAGCAACAATTTCTTACTTATGATAAAGGCACAATCTCCCC");

        final RunSGAViaProcessBuilderOnSpark.ContigsCollection contigsCollection = new RunSGAViaProcessBuilderOnSpark.ContigsCollection(contigsData);

        final List<ContigAligner.AssembledBreakpoint> assembledBreakpoints = contigAligner.alignContigs(contigsCollection);
        Assert.assertEquals(2, assembledBreakpoints.size());

        final AssembledBreakpoint breakpoint1 = assembledBreakpoints.get(0);

        Assert.assertEquals(breakpoint1.contigId, ">contig-9 2009 0");

        final AlignmentRegion breakpoint1Region1 = breakpoint1.region1;
        //todo verify that the 1-based vs 0-based nature of this is correct
        Assert.assertEquals(breakpoint1Region1.referenceInterval, new SimpleInterval("21", 27373209, 27374159));
        Assert.assertTrue(breakpoint1Region1.forwardStrand);
        Assert.assertEquals(breakpoint1Region1.mqual, 60);

        final AlignmentRegion breakpoint1Region2 = breakpoint1.region2;
        Assert.assertEquals(breakpoint1Region2.referenceInterval, new SimpleInterval("21", 27374159, 27374707));
        Assert.assertFalse(breakpoint1Region2.forwardStrand);
        Assert.assertEquals(breakpoint1Region2.mqual, 60);

        Assert.assertEquals(breakpoint1.homology, "GGATCCA");
        Assert.assertEquals(breakpoint1.insertedSequence, "NA");

        final AssembledBreakpoint breakpoint2 = assembledBreakpoints.get(1);

        Assert.assertEquals(breakpoint2.contigId, ">contig-9 2009 0");

        final AlignmentRegion breakpoint2Region1 = breakpoint2.region1;
        Assert.assertEquals(breakpoint2Region1.referenceInterval, new SimpleInterval("21", 27374159, 27374707));
        Assert.assertFalse(breakpoint2Region1.forwardStrand);
        Assert.assertEquals(breakpoint2Region1.mqual, 60);

        final AlignmentRegion breakpoint2Region2 = breakpoint2.region2;
        Assert.assertEquals(breakpoint2Region2.referenceInterval, new SimpleInterval("21", 27374701, 27375219));
        Assert.assertTrue(breakpoint2Region2.forwardStrand);
        Assert.assertEquals(breakpoint2Region2.mqual, 60);

        Assert.assertEquals(breakpoint2.homology, "NA");
        Assert.assertEquals(breakpoint2.insertedSequence, "NA");


    }

    @Test
    public void testAlignContigs2() throws Exception {

        // data taken from G94982 NA12878 assembly for breakpoint interval 14591
        final List<String> contigsData = new ArrayList<>();
        contigsData.add(">contig-3 312 0");
        contigsData.add("CCTGTAGATAGAGAGGTGGGTGAGAGATGGCCTTGTGGCAGCTCCTGGCAAGCTCACCTGACTTCTCATGATCTGGGTGACCATGGGGTATCCCTCCAAGACTTAGGTCAGCAGTGGTTAAGCCTTGCCCTGTAGCCTAGGAAAAAATGTGCAAGGTTGTCAGGGCACCAGCATGGAGGAGTTCCCCTACAGTCTTTCCAATACCTATGTGGTCTCTGGAACAGACATTTCATCCAGTAGCCATTCCTTTCCATTGTTTCCCTTCTTGGAAGAGCCTATCTTCCAAGACAGATGGTGAAATATTAGTAATTT");

        final RunSGAViaProcessBuilderOnSpark.ContigsCollection contigsCollection = new RunSGAViaProcessBuilderOnSpark.ContigsCollection(contigsData);

        final List<ContigAligner.AssembledBreakpoint> assembledBreakpoints = contigAligner.alignContigs(contigsCollection);
        Assert.assertEquals(1, assembledBreakpoints.size());

        final AssembledBreakpoint breakpoint1 = assembledBreakpoints.get(0);

        Assert.assertEquals(breakpoint1.contigId, ">contig-3 312 0");

        final AlignmentRegion breakpoint1Region1 = breakpoint1.region1;
        Assert.assertEquals(breakpoint1Region1.referenceInterval, new SimpleInterval("20", 1388956, 1389147));
        Assert.assertTrue(breakpoint1Region1.forwardStrand);
        Assert.assertEquals(breakpoint1Region1.mqual, 60);
        Assert.assertEquals(breakpoint1Region1.startInAssembledContig, 1);
        Assert.assertEquals(breakpoint1Region1.endInAssembledContig, 191);

        final AlignmentRegion breakpoint1Region2 = breakpoint1.region2;
        Assert.assertEquals(breakpoint1Region2.referenceInterval, new SimpleInterval("20", 1390815, 1390939));
        Assert.assertTrue(breakpoint1Region2.forwardStrand);
        Assert.assertEquals(breakpoint1Region2.mqual, 60);
        Assert.assertEquals(breakpoint1Region2.startInAssembledContig, 189);
        Assert.assertEquals(breakpoint1Region2.endInAssembledContig, 312);

        Assert.assertEquals(breakpoint1.homology, "ACA");
        Assert.assertEquals(breakpoint1.insertedSequence, "NA");
    }

    @Test
    public void testAlignedBreakpointBreakpointAlelle() throws Exception {
        final AlignmentRegion region1 = new AlignmentRegion(TextCigarCodec.decode("146M51S"), true, new SimpleInterval("8", 108569148, 108569294), 60, 1, 146, 0);
        final AlignmentRegion region2 = new AlignmentRegion(TextCigarCodec.decode("147S50M"), false, new SimpleInterval("8", 108569314, 108569364), 60, 148, 197, 0);
        final AssembledBreakpoint assembledBreakpoint = new AssembledBreakpoint("contig-1", region1, region2, "TC", "", new ArrayList<>());
        final SimpleInterval leftAlignedLeftBreakpointOnAssembledContig = assembledBreakpoint.getLeftAlignedLeftBreakpointOnAssembledContig();
        Assert.assertEquals(leftAlignedLeftBreakpointOnAssembledContig, new SimpleInterval("8", 108569294, 108569294));
        final SimpleInterval leftAlignedRightBreakpointOnAssembledContig = assembledBreakpoint.getLeftAlignedRightBreakpointOnAssembledContig();
        Assert.assertEquals(leftAlignedRightBreakpointOnAssembledContig, new SimpleInterval("8", 108569364, 108569364));
    }

    @Test
    public void testGetAssembledBreakpointsFromAlignmentRegions() throws Exception {
        final byte[] contigSequence = "GACGAACGATTTGACTTTAATATGAAATGTTTTATGTGGGCTATAAAATTATCCAAACTCGACACAGGACATTTTGAGCTTATTTCCAAATCATCTGGCCTTCATCTACCCACTGGAACTATTACTCTGCTGGGTCCTCATGGAAACATATCTTTCAGCCCTAACAATGAGACTACAGACATCTACGTCCCCAACACAACAGCTAAAAAGCAGTAGAATGTCAGAAAGGCTATCCACTTAGCCCTTGGCTGACAGGCCCCACTGAGCATCCTTTGCGAAGTCCATTTACTAGCTAATTCATAATTTACACAAGGCATTCAGACATAGCAGCTAAGATATAAAACATTTATCAACACAGGGACTAGTTTGTCATTTTAAAATAATTATGTTTAAGTAAGCCAATAAAGTCTATCTTCTCCAATTTACTTATTGAGCTTTATGAGGCAATTTAAGTCCCGATTTTGGGGGGTATGTATGAAAGGAGAGCATGGAAATGCCATTTGCTCCCTGAAGTTTTTATCTTTTTTTTTTTGAGATAGAGTCTTGTGTTTTCTGTGGAGTACATGAGTATGCATCAAAGCTAACAACGCCCACTGCCCTGTTAGTCAAATACCTTTGA".getBytes();
        final AlignmentRegion region1 = new AlignmentRegion(TextCigarCodec.decode("532M87S"), true, new SimpleInterval("8", 118873207, 118873739), 60, 1, 532, 0);
        final AlignmentRegion region2 = new AlignmentRegion(TextCigarCodec.decode("518S29M72S"), false, new SimpleInterval("1", 175705642, 175705671), 3, 519, 547, 0);
        final AlignmentRegion region3 = new AlignmentRegion(TextCigarCodec.decode("543S76M"), false, new SimpleInterval("1", 118875262, 118875338), 60, 544, 619, 0);
        final ArrayList<AlignmentRegion> alignmentRegionList = new ArrayList<>();
        alignmentRegionList.add(region1);
        alignmentRegionList.add(region2);
        alignmentRegionList.add(region3);
        final List<AssembledBreakpoint> assembledBreakpointsFromAlignmentRegions = ContigAligner.getAssembledBreakpointsFromAlignmentRegions("contig-0", contigSequence, alignmentRegionList);
        Assert.assertEquals(assembledBreakpointsFromAlignmentRegions.size(), 1);
        final AssembledBreakpoint assembledBreakpoint = assembledBreakpointsFromAlignmentRegions.get(0);
        Assert.assertEquals(assembledBreakpoint.contigId, "contig-0");
        Assert.assertEquals(assembledBreakpoint.region1, region1);
        Assert.assertEquals(assembledBreakpoint.region2, region3);
        Assert.assertEquals(assembledBreakpoint.homology, "NA");
        Assert.assertEquals(assembledBreakpoint.insertedSequence, "GAGATAGAGTCT");
    }

    @Test
    public void testGetAssembledBreakpointsFromAlignmentRegionsWithOverlappingAlignmentRegion() throws Exception {
        final byte[] contigSequence = "ACTAGAGCATCTACGTGTTCCTGTGGTTTTGGAGCAAGAGTGATTTGAGTTTCAGAGATTTTTACTAATTCTTCTTCCCCTACCAGAAAAAAAGATCTTACCATTTGAGAGTGAGATGTAAACCCAGCCCTGTCTGACCTGAGTCTGTGCCCTAAGCCTATGCTAAGCCAAGCAGTGCCTGGAGCCACCACAGGTCCACACAATTCGTTAACATGATGAAGCAAGGATGGAAATTGGACAAAATAGTGTGCCTACTGAATCTAAGAATGAAAAATGATTGCACTCCTACTCTGAGTGCTTTGGAGCACTGCCCAGTTGGGCAAAGGGTCAGCGCCTGGGCAGAGGTCCCCACAACCTGGCAGGAGTGTGGTCGGCCACCCTATGGGCCTCCATCATGTGCAGTGACAGCGGGGCTGTCATGTCACCGTGTGGGAGGGCTTGCAGGTGAAGTGGTCTGGGAGGGGTCCCCCAGACAAAGCCAAGGTTCTGAGAGTTGGCCCGAACACTGCTGGATTCCACTTCACCTGCAAGCCCTCCCACACGGTGACATGACAGCCTATAATACAGTTCCGCATGGCCACGTCATACAACCCTGTCATATTGGTGAGCAATTGCTGTGTAGCCAAAGACCCCAAAACTCAAACAGCATTTATTATTATTGCCCCCATGTCTGAGAGTCAGATGTGCATTTGCTGATCTCAGCTTGTTTGAGCTGCTGCAGGGTTGGGGCTCTGCTCCAGGCAGGCTTAGCTGTCACCACATGCACACATACATTCTGGGCCTCTGCTGCGCGCGTCACGTTCACTGAAGATCTTGGGATTGGGAGTTAGGGCGGTGGGAGGGCCCAGCAAAGTCACCTGGCGATGGCAGGGACACAGGGAGGAATGTAGAATGGGGCCGATGATGGGACCCACACGTCTGCAAAGCTGCGGTCTCCTTGAGGGGTGGAGACAGCAACAACTCACCGCACGCGGTGCTTCAGTTCACCATCTCCCTGGGACATTAGGGGGCCCCGTGTTATCTCATTTTGCTCTGGTTTGCATTAGTTTTTTATCACTTCGTAGATGAAGCCACTGACACCCAGAGAGGGAAAGTGGCCTGACCAAGGGCCACAGCAGGGGAGCGAAGGAGCCCCACAGTTCGGCAGGAACACAGCCTCTCCCTGGCTTTCAGGTTCACTGACATCTTCTCATGGCCTCTGTAACTCACCAGGCATCAGGGTGTAGTCCTTAGACCAGTGTCCCACAGCTGCCACAGAGTGGGAGCTCACCATCAGTTATAAGTCACTAGAAAGGCTTTTGGACATTATAAGCTACAATGGAAAATAAGTCATCTGTGGATTTTTGTGACAGATTCCAAAAATTTGAATATTTTGTCTACTTAGGTTTTTGGTTAATTTTATCCTCAAAACTGTTCTGCAGTGATTAAGCTGTACAAACTGCATCATGGGCGAATTGGCATATTCAGAAATGACTGATATTCTTGATTTCAGTTTTTTACTTTGTATGTAGCTCCTCAAGGAAAC".getBytes();
        final AlignmentRegion region1 = new AlignmentRegion(TextCigarCodec.decode("487M1006S"), true, new SimpleInterval("20", 23102817, 23103304), 60, 1, 487, 1);
        final AlignmentRegion region2 = new AlignmentRegion(TextCigarCodec.decode("483S42M968S"), false, new SimpleInterval("20", 23103196, 23103238), 60, 484, 525, 2);
        final AlignmentRegion region3 = new AlignmentRegion(TextCigarCodec.decode("523S970M"), true, new SimpleInterval("20", 23103633, 23104603), 60, 524, 1493, 3);
        final ArrayList<AlignmentRegion> alignmentRegionList = new ArrayList<>();
        alignmentRegionList.add(region1);
        alignmentRegionList.add(region2);
        alignmentRegionList.add(region3);
        final List<AssembledBreakpoint> assembledBreakpointsFromAlignmentRegions = ContigAligner.getAssembledBreakpointsFromAlignmentRegions("contig-0", contigSequence, alignmentRegionList);
        Assert.assertEquals(assembledBreakpointsFromAlignmentRegions.size(), 1);
        final AssembledBreakpoint assembledBreakpoint = assembledBreakpointsFromAlignmentRegions.get(0);
        Assert.assertEquals(assembledBreakpoint.contigId, "contig-0");
        Assert.assertEquals(assembledBreakpoint.region1, region1);
        Assert.assertEquals(assembledBreakpoint.region2, region3);
        Assert.assertEquals(assembledBreakpoint.homology, "NA");
        Assert.assertEquals(assembledBreakpoint.insertedSequence, "TGAGAGTTGGCCCGAACACTGCTGGATTCCACTTCAC");
        Assert.assertEquals(assembledBreakpoint.insertionMappings.size(), 1);
        Assert.assertEquals(assembledBreakpoint.insertionMappings.get(0), "484-525:20,23103196,-,483S42M968S,60,2");
    }


    @AfterClass
    public void tearDown() throws Exception {
        contigAligner.close();
    }
}
