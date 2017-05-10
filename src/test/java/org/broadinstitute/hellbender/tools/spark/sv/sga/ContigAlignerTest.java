package org.broadinstitute.hellbender.tools.spark.sv.sga;

import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedAssembly;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexSingleton;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class ContigAlignerTest extends BaseTest {

    private ContigAligner contigAligner;

    @BeforeClass
    public void setup() throws Exception {
        contigAligner = new ContigAligner(b37_reference_20_21+".img");
    }

    @Test(groups = "sv")
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

        final ContigsCollection contigsCollection = new ContigsCollection(contigsData);

        final AlignedAssembly alignedAssembly = contigAligner.alignContigs(1, contigsCollection);
        Assert.assertEquals(alignedAssembly.alignedContigs.size(), 18);
        Assert.assertEquals(alignedAssembly.alignedContigs.stream().mapToInt(ctg -> ctg.alignmentIntervals.size()).sum(), 20);

        final AlignedContig alignedContig_9 = alignedAssembly.alignedContigs.get(3);
        final List<AlignedAssembly.AlignmentInterval> alignmentIntervalsForContig9 = alignedContig_9.alignmentIntervals;
        Assert.assertEquals(alignedContig_9.contigName, "asm000001:tig00009");
        Assert.assertEquals(alignmentIntervalsForContig9.size(), 3);

        final AlignedAssembly.AlignmentInterval firstAlignmentIntervalForContig_9 = alignmentIntervalsForContig9.get(0);
        Assert.assertEquals(firstAlignmentIntervalForContig_9.referenceInterval, new SimpleInterval("21", 27373209, 27374158));
        Assert.assertTrue(firstAlignmentIntervalForContig_9.forwardStrand);
        Assert.assertEquals(firstAlignmentIntervalForContig_9.mapQual, 60);

        final AlignedAssembly.AlignmentInterval secondAlignmentIntervalForContig_9 = alignmentIntervalsForContig9.get(1);
        Assert.assertEquals(secondAlignmentIntervalForContig_9.referenceInterval, new SimpleInterval("21", 27374159, 27374706));
        Assert.assertFalse(secondAlignmentIntervalForContig_9.forwardStrand);
        Assert.assertEquals(secondAlignmentIntervalForContig_9.mapQual, 60);

        final AlignedAssembly.AlignmentInterval thirdAlignmentIntervalForContig_9 = alignmentIntervalsForContig9.get(2);
        Assert.assertEquals(thirdAlignmentIntervalForContig_9.referenceInterval, new SimpleInterval("21", 27374701, 27375218));
        Assert.assertTrue(thirdAlignmentIntervalForContig_9.forwardStrand);
        Assert.assertEquals(thirdAlignmentIntervalForContig_9.mapQual, 60);
    }

    @Test(groups = "sv")
    public void testAlignContigs2() throws Exception {

        // data taken from G94982 NA12878 assembly for breakpoint interval 14591
        final List<String> contigsData = new ArrayList<>();
        contigsData.add(">contig-3 312 0");
        contigsData.add("CCTGTAGATAGAGAGGTGGGTGAGAGATGGCCTTGTGGCAGCTCCTGGCAAGCTCACCTGACTTCTCATGATCTGGGTGACCATGGGGTATCCCTCCAAGACTTAGGTCAGCAGTGGTTAAGCCTTGCCCTGTAGCCTAGGAAAAAATGTGCAAGGTTGTCAGGGCACCAGCATGGAGGAGTTCCCCTACAGTCTTTCCAATACCTATGTGGTCTCTGGAACAGACATTTCATCCAGTAGCCATTCCTTTCCATTGTTTCCCTTCTTGGAAGAGCCTATCTTCCAAGACAGATGGTGAAATATTAGTAATTT");

        final ContigsCollection contigsCollection = new ContigsCollection(contigsData);

        final AlignedAssembly alignedAssembly = contigAligner.alignContigs(1, contigsCollection);
        Assert.assertEquals(alignedAssembly.alignedContigs.size(), 1);
        Assert.assertEquals(alignedAssembly.alignedContigs.stream().mapToInt(ctg -> ctg.alignmentIntervals.size()).sum(), 2);

        final AlignedContig alignedContig_3 = alignedAssembly.alignedContigs.get(0);
        final List<AlignedAssembly.AlignmentInterval> alignmentIntervalsForContig_3 = alignedContig_3.alignmentIntervals;
        Assert.assertEquals(alignedContig_3.contigName, "asm000001:tig00003");
        Assert.assertEquals(alignmentIntervalsForContig_3.size(), 2);

        final AlignedAssembly.AlignmentInterval firstAlignmentIntervalForContig_3 = alignmentIntervalsForContig_3.get(0);
        Assert.assertEquals(firstAlignmentIntervalForContig_3.referenceInterval, new SimpleInterval("20", 1388956, 1389146));
        Assert.assertTrue(firstAlignmentIntervalForContig_3.forwardStrand);
        Assert.assertEquals(firstAlignmentIntervalForContig_3.mapQual, 60);
        Assert.assertEquals(firstAlignmentIntervalForContig_3.startInAssembledContig, 1);
        Assert.assertEquals(firstAlignmentIntervalForContig_3.endInAssembledContig, 191);

        final AlignedAssembly.AlignmentInterval secondAlignmentIntervalForContig_3 = alignmentIntervalsForContig_3.get(1);
        Assert.assertEquals(secondAlignmentIntervalForContig_3.referenceInterval, new SimpleInterval("20", 1390815, 1390938));
        Assert.assertTrue(secondAlignmentIntervalForContig_3.forwardStrand);
        Assert.assertEquals(secondAlignmentIntervalForContig_3.mapQual, 60);
        Assert.assertEquals(secondAlignmentIntervalForContig_3.startInAssembledContig, 189);
        Assert.assertEquals(secondAlignmentIntervalForContig_3.endInAssembledContig, 312);
    }

    @Test(groups = "sv")
    public void testAlignArtificialBreakpointContig() throws Exception {

        // take DNA sequence from the reference split between 20:
        final List<String> contigsData = new ArrayList<>();
        contigsData.add(">contig-20 fake 1000000-1000099+20:5000002-5000101");// old format that will throw exception: contigsData.add(">contig-fake-20:1000000-1000099+20:5000002-5000101");
        final String snippet1 = "GTGGGAGAGAACTGGAACAAGAACCCAGTGCTCTTTCTGCTCTACCCACTGACCCATCCTCTCACGCATCATACACCCATACTCCCATCCACCCACCTTC";
        final String snippet2 = "GTGATCCAGCTACAGACTGTTCCAAAGACTTTGCAACTGTTATTTTTGCTTAATCCTCACAACAACCTATGAGGTAGGCACATTTATTGCCCCCATGTGA";
        contigsData.add(snippet1 + snippet2);

        final ContigsCollection contigsCollection = new ContigsCollection(contigsData);

        final AlignedAssembly alignedAssembly = contigAligner.alignContigs(21, contigsCollection);
        Assert.assertEquals(alignedAssembly.alignedContigs.size(), 1);
        Assert.assertEquals(alignedAssembly.alignedContigs.stream().mapToInt(ctg -> ctg.alignmentIntervals.size()).sum(), 2);

        final AlignedContig alignedContig_20 = alignedAssembly.alignedContigs.get(0);
        final List<AlignedAssembly.AlignmentInterval> alignmentIntervalsForContig_20 = alignedContig_20.alignmentIntervals;
        Assert.assertEquals(alignedContig_20.contigName, "asm000021:tig00020");
        Assert.assertEquals(alignmentIntervalsForContig_20.size(), 2);

        final AlignedAssembly.AlignmentInterval firstAlignmentIntervalForContig_20 = alignmentIntervalsForContig_20.get(0);
        Assert.assertEquals(firstAlignmentIntervalForContig_20.referenceInterval, new SimpleInterval("20", 1000000, 1000099));
        Assert.assertTrue(firstAlignmentIntervalForContig_20.forwardStrand);
        Assert.assertEquals(firstAlignmentIntervalForContig_20.mapQual, 60);
        Assert.assertEquals(firstAlignmentIntervalForContig_20.startInAssembledContig, 1);
        Assert.assertEquals(firstAlignmentIntervalForContig_20.endInAssembledContig, 100);

        final AlignedAssembly.AlignmentInterval secondAlignmentIntervalForContig_20 = alignmentIntervalsForContig_20.get(1);
        Assert.assertEquals(secondAlignmentIntervalForContig_20.referenceInterval, new SimpleInterval("20", 5000002, 5000101));
        Assert.assertTrue(secondAlignmentIntervalForContig_20.forwardStrand);
        Assert.assertEquals(secondAlignmentIntervalForContig_20.mapQual, 60);
        Assert.assertEquals(secondAlignmentIntervalForContig_20.startInAssembledContig, 101);
        Assert.assertEquals(secondAlignmentIntervalForContig_20.endInAssembledContig, 200);
    }

    @AfterClass
    public void tearDown() throws Exception {
        BwaMemIndexSingleton.closeInstance();
    }
}
