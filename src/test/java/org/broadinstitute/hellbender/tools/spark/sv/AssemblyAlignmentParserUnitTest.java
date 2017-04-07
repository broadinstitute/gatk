package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.collect.Iterators;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.apache.commons.collections4.IterableUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;

import static org.testng.Assert.assertEquals;


public class AssemblyAlignmentParserUnitTest extends BaseTest{
    private static final String dummyRefName = "1";
    private static final int dummyRefId = Integer.valueOf(dummyRefName) - 1;
    private static final List<String> refNames = Collections.singletonList(dummyRefName);


    @Test
    public void testConvertGATKReadToAlignmentRegions() throws Exception {
        final String read1Seq = "ACACACACACACACACACACACACACACCAGAAGAAAAATTCTGGTAAAACTTATTTGTTCTTAAACATAAAACTAGAGGTGCAAAATAACATTAGTGTATGGTTAATATAGGAAAGATAAGCAATTTCATCTTCTTTGTACCCCACTTATTGGTACTCAACAAGTTTTGAATAAATTCGTGAAATAAAGAAAAAAACCAACAAGTTCATATCTCCAAGTAGTCTTGGTAATTTAAACACTTTCAAGTATCTTCATACCCTGTTAGTACTTCTTTCTGCTCTGTGTCAACAGAAGTATTCTCAACACTCTTGTGGTTAATTTGGTTAAAACTCATTACAAAGAACCCTAAGTTATCCGTCACAGCTGCTAAACCCATTAGTACAGGAGACTTTAAGGCCATAATGTGTCCATTTTCAGGATATAATTGAAGAGAGGCAAATGATACATGGTTTTCCAAAAATATTGGACCAGGGAGCCTCTTCAAGAAAGAATCCCTGATTCGGGAGTTCTTATACTTTTTCAAGAA";
        final byte[] read1Quals = new byte[read1Seq.length()];
        Arrays.fill(read1Quals, (byte)'A');
        final String read1Cigar = "527M2755H";

        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(read1Seq.getBytes(),
                read1Quals, read1Cigar);
        read1.setPosition("chrX", 6218482);
        read1.setIsReverseStrand(true);
        read1.setIsSupplementaryAlignment(true);

        final String read2Seq = "ATATATATATATATATATATATATATAATATAAGAAGGAAAAATGTGGTTTTCCATTATTTTCTTTTTCTTATCCTCATCATTCACCAAATCTATATTAAACAACTCATAACATCTGGCCTGGGTAATAGAGTGAGACCCCAACTCCACAAAGAAACAAAAATTAAAAACAAATTAGCCGGGCCTGATGGCAAGTACCTGTGGTTCCAGCTGAGGTGGGAGGATCACTTGAGCCCAGGAGTTCAGGGCTGCAGTGAGCTATGATTACGCCAGTGTACTCCAGCCTGGGAGACAGAGCAAGACCCTATCTCTAAAAATATAAATACATAAATAAATAAATAATAAATAAAAATAGAAAATGTACAATGAAAGTTATAAAGTTGGCCAGGCGTGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCTGAGGTGGGCGAATCACCTGAGGTCAGTAGTTCAAGACTAGCCTGGCCAACATGGCGAAATCCTGTCTCTACTAAAAATACAAAAACTAGCTGGGTGTGGTGGTGTGTGCCTGTAATCCCAGCTATACAGGAGGCTGAGGCCGGAGAATTGCTTGAACCTGGGAGGGGGAGGTTGCAGTGAGCCAAGATCGTGCCATTGCACTGTAGCTTGGGTGACAGAGCGAGACTCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAGTTATAAAGTTACCTATGATGGGTCTGGATGTACTCCTTATTTAGGAGTGAAGACATTCGTTAACATGAGACCTAAGTAAGTAGAAAGTATGTGTTTAAGGGACAGGTGTCCATTTTCTCTAGGTCTCCTGGAAGCTTTTTTTTTCTAATTTGAGTACTAGTTCCAAAAAAGGTGTTACCGCCTATGTTTATAGTGAAACTATCTATGTGTGACAAAATTCTACCCTCTCTTGTCCATCAATATTGTGCAATGTTGTGTACTTGTATGGAGAAATAGACAACTTTTACAAAGATCAAACTAGGCACCCTTTACCAACGCTAAACTCATAACCCTTTTATCTGCCTTTGTAGAAGATTCTCACCTTTATTTCTCTTGGTCCCCAAAGGCTTCGCTGAGAGACTGTTCCCATTGAAGCCTTGTGGCAAAGTCAGTAGAGCATTGTATCAGCCCAGCTTCAGAAACTCTGCATTCACTTTTAATAAATATGGAGAAGTTTGAAGTCACAAAAGCTGAGAACCTCATTACAGTCTCTTTATGTTCTTGGGTCCCTCTCTTCCTGCAGACTTCCTTGAATTCAAATTTAATGTGCGCTACTCATTCACATTGCACCTACTGAAAAGCATTGACAATTCCAGGTGACGACAGGAGCAAAAGGGAATGAATGCTAGAGGGTTTCTCATTTAGGTTGTTCACTCACCACTGTGCACAGGTTCTTTAGAATTCTACCTACACATGAAATAATATATAGTAAAATCCATATACATCTATAATAAATATATTCCAATGTATATATATTGTATATATATAGTTCTTGTTCTTTTGTATGAAAAGTACATGCAAATTTTATATGTATATACATATATACGTATATATATGTATGTGTGTGTGTGTATATATATATATACACACACACAAACCATATATATAGTTTGACATGCATTTTCCAAAAGAAACTATATTGTATGTGGTGAATGAACTTCCAGGCTGGTTCACAGGAGCTATTGTTAGTTAATGGTGCACAGAAACCACAGCTTCCTTCCCTCTTCCCTTCCTTTCTTCTTTCCTTCCTTCCTTCCCTCCTTCCTTCCTTCTTTTTCTCTTCTTTCTTTCCTTCCTTCTTTCCTTCATTCCCTCCCTCTCTCCTTCCTACCTTTTTCTCTTCTTTCTTTCCTTCCTTTCTTCTTTCCCTCCTTCCCTCCTTCTCTCCTTCCCTTCGATCTATCATCCTTCTCCCTCCCTCCCTGCCTCCTTTCTTCCTTCCTTTTCTCTTTTCTTTCTTTCTTCCTTCCTTCCTTCCTCCTTTCCTTCCCTTCCTTCCTCCCTTCGTTCCTCCATCCCTTCCTTCTTCCTTCTGTCCTTTCTTCCTTCCTTCTTCCCTTCCTTCTTTCTTCCCTCATTCCTTTATTCTATTTGGATGAGTGACAAGGTCTCAGGGATCCTTAAAGTATTAAAGAGAAAATGAAGATAATTCTATAGAATAACAGATCAAATAGGAAAGAGTAATGAAAAGAAAGACTCCTGGAATGACAAAGAAGCATGTTCATGACAATTTAATAGAATTTCCCACATGCCATATTTGTATAGCTGGTTATATTTTAAAACTTTGTCAAACCCTATATTTATAAAACCTTGATTTTGATTTTAGTAAATTCTATGAATAAATGCATACATAAAGCAGATTACTTTGATTCGACAAAAACATACAGGCACCTGTGAGCTTTTTTTTTTTTTATGTGGTAGTGGAAACAGCCTCATTAATTAATGGCTATTTCCTATATTCTTTTTAAAATTAATAAGCAACAGTCAGATATTTCAGCATTTTTAAGGCATGAACAATACCTTGGATTCTAATATATTTTAGTCTAAAATATCTGGTAACATCACAGATAATTGTTTCAAACTCTACCACCCCTCACACAATTTATTTTTCAACAATTTAAAGTCTTTTAGAATAATGAGATAACTTTAAAAATAGACAACATTACTCATACTGGTGAGCTAGGATTGCTTAAAGATGGGGCATAGATTGCATTGTCTCAGAGGAAAATATTTTATGAAGATTCTTGAAAAAGTATAAGAACTCCCGAATCAGGGATTCTTTCTTGAAGAGGCTCCCTGGTCCAATATTTTTGGAAAACCATGTATCATTTGCCTCTCTTCAATTATATCCTGAAAATGGACACATTATGGCCTTAAAGTCTCCTGTACTAATGGGTTTAGCAGCTGTGACGGATAACTTAGGGTTCTTTGTAATGAGTTTTAACCAAATTAACCACAAGAGTGTTGAGAATACTTCTGTTGACACAGAGCAGAAAGAAGTACTAACAGGGTATGAAGATACTTGAAAGTGTTTAAATTACCAAGACTACTTGGAGATATGAACTTGTTGGTTTTTTTCTTTATTTCACGAATTTATTCAAAACTTGTTGAGTACCAATAAGTGGGGTACAAAGAAGATGAAATTGCTTATCTTTCCTATATTAACCATACACTAATGTTATTTTGCACCTCTAGTTTTATGTTTAAGAACAAATAAGTTTTACCAGAATTTTTCTTCTGGTGTGTGTGTGTGTGTGTGTGTGTGTGT";
        final byte[] read2Quals = new byte[read2Seq.length()];
        Arrays.fill(read2Quals, (byte)'A');
        final String read2Cigar = "1429S377M8D34M4I120M16D783M535S";

        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(read2Seq.getBytes(),
                read2Quals, read2Cigar);
        read2.setPosition("chrX", 6219006);

        final String read3Seq = "GGGACCAAGAGAAATAAAGGTGAGAATCTTCTACAAAGGCAGATAAAAGGGTTATGAGTTTAGCGTTGGTAAAGGGTGCCTAGTTTGATCTTTGTAAAAGTTGTCTATTTCTCCATACAAGTACACAACATTGCACAATATTGATGGACAAGAGAGGGTAGAATTTTGTCACACATAGATAGTTTCACTATAAACATAGGCGGTAACACCTTTTTTGGAACTAGTACTCAAATTAGAAAAAAAAAGCTTCCAGGAGACCTAGAGAAAATGGACACCTGTCCCTTAAACACATACTTTCTACTTACTTAGGTCTCATGTTAACGAATGTCTTCACTCCTAAATAAGGAGTACATCCAGACCCATCATAGGTAACTTTATAACTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACAGAGTCTCGCTCTGTCACCCAAGCTACAGTGCAATGGCACGATCTTGGCTCACTGCAACCTCCCCCTCCCAGGTTCAAGCAATTCTCCGGCCTCAGCCTCCTGTATAGCTGGGATTACAGGCACACACCACCACACCCAGCTAGTTTTTGTATTTTTAGTAGAGACAGGATTTCGCCATGTTGGCCAGGCTAGTCTTGAACTACTGACCTCAGGTGATTCGCCCACCTCAGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCACGCCTGGCCAACTTTATAACTTTCATTGTACATTTTCTATTTTTATTTATTATTTATTTATTTATGTATTTATATTTTTAGAGATAGGGTCTTGCTCTGTCTCCCAGGCTGGAGTACACTGGCGTAATCATAGCTCACTGCAGCCCTGAACTCCTGGGCTCAAGTGATCCTCCCACCTCAGCTGGAACCACAGGTACTTGCCATCAGGCCCGGCTAATTTGTTTTTAATTTTTGTTTCTTTGTGGAGTTGGGGTCTCACTCTATTACCCAGGCCAGATGTTATGAGTTGTTTAATATAGATTTGGTGAATGATGAGGATAAGAAAAAGAAAATAATGGAAAACCACATTTTTCCTTCTTATATTATATATATATATATATATATATATAT";
        final byte[] read3Quals = new byte[read3Seq.length()];
        Arrays.fill(read3Quals, (byte)'A');
        final String read3Cigar = "2207H385M6I684M";

        final GATKRead read3 = ArtificialReadUtils.createArtificialRead(read3Seq.getBytes(),
                read3Quals, read3Cigar);
        read3.setIsReverseStrand(true);
        read3.setIsSupplementaryAlignment(true);

        List<GATKRead> reads = new ArrayList<>(3);
        reads.add(read1);
        reads.add(read2);
        reads.add(read3);

        final Tuple2<Iterable<AlignmentRegion>, byte[]> alignmentRegionResults = AssemblyAlignmentParser.convertToAlignmentRegions(reads);
        assertEquals(alignmentRegionResults._2(), read2.getBases());

        final List<AlignmentRegion> alignmentRegions = IterableUtils.toList(alignmentRegionResults._1);
        assertEquals(alignmentRegions.size(), 3);

        final byte[] read4Bytes = SVCallerTestDataProvider.LONG_CONTIG1.getBytes();
        final String read4Seq = new String(read4Bytes);
        final byte[] read4Quals = new byte[read4Seq.length()];
        Arrays.fill(read4Quals, (byte)'A');
        final String read4Cigar = "1986S236M2D1572M1I798M5D730M1I347M4I535M";

        final GATKRead read4 = ArtificialReadUtils.createArtificialRead(read4Seq.getBytes(),
                read4Quals, read4Cigar);
        read4.setIsReverseStrand(true);
        read4.setIsSupplementaryAlignment(false);

        final String read5Seq = "TTTCTTTTTTCTTTTTTTTTTTTAGTTGATCATTCTTGGGTGTTTCTCGCAGAGGGGGATTTGGCAGGGTCATAGGACAACAGTGGAGGGAAGGTCAGCAGATAAACAAGTGAACAAAGGTCTCTGGTTTTCCTAGGCAGAGGACCCTGCGGCCTTCCGCAGTGTTTGTGTCCCCGGGTACTTGAGATTAGGGAGTGGTGATGACTCTTAACGAGCATGCTGCCTTCAAGCATCTGTTTAACAAAGCACATCTTGCACCGCCCTTAATCCATTTAACCCTGAGTGGACACAGCACATGTTTCAGAGAGCACAGGGTTGGGGGTAAGGTCATAGGTCAACAGGATCCCAAGGCAGAAGAATTTTTCTTAGTACAGAACAAAATGAAGTCTCCCATGTCTACCTCTTTCTACACAGACACAGCAACCATCCGATTTCTCAATCTTTTCCCCACCTTTCCCCCGTTTCTATTCCACAAAACTGCCATTGTCATCATGGCCCGTTCTCAATGAGCTGTTGGGTACACCTCCCAGACGGGGTGGTGGCCGGGCAGAGGGGCTCCTCACTTCCCAGTAGGGGCGGCCGGGCAGAGGCACCCCCCACCTCCCGGACGGGGCGGCTGGCTGGGCGGGGGGCTGACCCCCCCCACCTCCCTCCCGGACGGGGCGTCAGCCCAGCGTTTCTGATTGGATAATGCTTAAGGCCCCCGCCCCCTCAGGCCCTGAGTAACAGAAAATGTGATCAGGACTGAGTGAAGAAAAAGTCACAGCCTAAGCTGCAGCGTTTTTCAGGCAGGGCTTCCTCCCTGAGCTAAGCCAGGCCCAACCCAGATCATGGGAAAATTCTATCTTTTTTTTTACTCTCTCTCTTTTTGAATATATTCAAAAGGTGAACAGAAGTATTTTGCTGTCATATTAATAATACATAAAATTTTTTTCAAGAGAAAATCAGCTTTTACTTTGGTAATAGTGTATTATCAATACTAAAGCTAATTTTAATAAACCTTATAAATAAATCAAATTTGTCATTTTTGACCACTCCGGTTTTACATGTATATTTTGTAATCTCTTGTAATTTTTAAAAACTGTTTACATTTTATTTTTATCCAAATTCTTTTTATTTTTTCAATTTGAAACCACCTTTAAGTAATTTCAAATTGTTATAGGAGATAGAAAGAAGTCATTTAGGGCCAGGTTCACTGGCAAGTGCCTGTAATCGCAACACTTTGGGAGGCCAAGGTGGGTGGATCACTTGAGGTCAGGAGTTCAAGACCAGCCTGGCCAACATGGTGAAACCCCATTTCTACTAAAAGTACAAATAACTAGCTGGGTGTGGTGCTGCACACCTGTAGACCCAGCTAGTCCGGAGGCTGAGGCAGGAGAATTGCTTAAACCCAGAAGGCGGAGGTTGTGTAACTGCCCAAGGGGTTCACCTTGCCCTATGTCTAGACAGAGCTGATTCATCAAGACTGGGTAATTTGTGGAGGAAAAGTTAAATACTAAATTTGAACTTAATTGAACGTGGACAAACTCAGTAGTCACCAAGTTCTCAAACAGGTTGTGTGAGGCCCTTGAGGCATTCATTCAGCGCTGTTTCAGAGAAATCTCTATTTCAATCTATTTCTATATATTAGTTGTTGAAAAACAATAGACAATCGCAAAAACAAGTTGACAGTTTTGTGTTCCTTGACCCCAGTTGCAAATGGCCCTCATGACTGGGCCTTATGCCAAACAACTCGTTACAAAAGAGCTAGGGTATCAGATTGTGCTGAAGCTTCATTAGACCCTCCTCATCTGTGCAGGAATGAGTAGCTGACTCTGGAGCCCAAGCTGTTGCTTCCCAGTCTGGTGGTGAATCCTCCATAGTCTGGTGAGTGTAAATATATATATCTCTTTTCCTTTCTCCTCTTCCCATTGCAATTTGCTTATTATATCAACATGCTTATTATATCAATCTGGTTATTACATTGATTTGCTTATTATATAATTTGCTAATTATATCTGCATTTCCATTTACGTGGCACAAAGCTTATTTACCCTTAAAGGTATTGTGTGTGTGTCTTTTCTTCTCCCCTTGAATGTTTCCCACACAGAATATTTTTGGCGTCACGAACAGGATTCAAAAACCAAACTGTGCCACTTTTTGGCCACAAGGACAGGGCTGGAAACTCGAGGAAATCCCAAATCCAGGGATGGGAACTCCCCAATATCACTGTCAATTCCTACAGGATTGAACGAAGGGGACGAATGCAGAAATGAAGACAAAGACAAAAGATTTGTTTTAAAAGAAGGGGTCAGGCAGGGCGCAGTGGCTCAGGCCTGTAATCCCAGCACTTTGAGAGGCCGAGGTGGGCGGATCGCGAGGTCAGGAGAGCGAAACCATCTTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAAATTTAGCCAGGTGTGGTGGCAGACACCTGTAGTCCCAGCTACCTGGGGGGGTGGGGGGGTGGGGCTGAGGCAGGAGAATGGCATGAACCCAGGAGGCAGAGCTTGCAGTAAGCCAAGATCGTGCCACTGCACTCCAGCCTGGGTGACAGAGCGAGATTCGGTCTCAGAAAAAAAAAAAAAAA";
        final byte[] read5Quals = new byte[read5Seq.length()];
        Arrays.fill(read5Quals, (byte)'A');
        final String read5Cigar = "3603H24M1I611M1I1970M";

        final GATKRead read5 = ArtificialReadUtils.createArtificialRead(read5Seq.getBytes(),
                read5Quals, read5Cigar);
        read5.setIsReverseStrand(false);
        read5.setIsSupplementaryAlignment(true);

        List<GATKRead> reads2 = new ArrayList<>(2);
        reads2.add(read4);
        reads2.add(read5);

        final Tuple2<Iterable<AlignmentRegion>, byte[]> alignmentRegionResults2 = AssemblyAlignmentParser.convertToAlignmentRegions(reads2);
        // these should be the reverse complements of each other
        assertEquals(alignmentRegionResults2._2().length, read4.getBases().length);

        final List<AlignmentRegion> alignmentRegions2 = IterableUtils.toList(alignmentRegionResults2._1);
        assertEquals(alignmentRegions2.size(), 2);

        final AlignmentRegion alignmentRegion4 = alignmentRegions2.get(0);
        assertEquals(alignmentRegion4.startInAssembledContig, 1);
        assertEquals(alignmentRegion4.endInAssembledContig, read4Seq.length() - 1986);

        final AlignmentRegion alignmentRegion5 = alignmentRegions2.get(1);
        assertEquals(alignmentRegion5.startInAssembledContig, 3604);
        assertEquals(alignmentRegion5.endInAssembledContig, read4Seq.length());
    }

    @Test
    public void testConvertAlignedAssemblyOrExcuseToAlignmentRegions() {

        // test "failed" assembly doesn't produce anything
        final AlignedAssemblyOrExcuse excuse = new AlignedAssemblyOrExcuse(1, "justATest");
        Assert.assertTrue(AssemblyAlignmentParser.formatToAlignmentRegions(Collections.singletonList(excuse), refNames).isEmpty());

        final byte[] dummyContigSequence = SVCallerTestDataProvider.makeDummySequence(1000, (byte)'T');
        final byte[] dummyContigSequenceQuals = SVCallerTestDataProvider.makeDummySequence(1000, (byte)'A');

        final FermiLiteAssembly.Contig unmappedContig = new FermiLiteAssembly.Contig(dummyContigSequence, dummyContigSequenceQuals, 100); // totally random 100 supporting reads
        final BwaMemAlignment unmappedContigAlignment = new BwaMemAlignment(4, -1, -1, -1, -1, -1, -1, -1, 0, 0, "", "", "", -1, -1, 0);

        final FermiLiteAssembly.Contig contigWithAmbiguousMapping = new FermiLiteAssembly.Contig(dummyContigSequence, dummyContigSequenceQuals, 100);
        final BwaMemAlignment firstAmbiguousMapping = new BwaMemAlignment(256, dummyRefId, 1000000, 1001000, 0, 1000, 0, 20, 100, 100, "800M50I100M50D50M", "", "", -1, -1, 0); // technically not correct but doesn't matter for this case
        final BwaMemAlignment secondAmbiguousMapping = new BwaMemAlignment(272, dummyRefId, 2000000, 2001000, 0, 1000, 0, 50, 100, 100, "700M50I200M50D50M", "", "", -1, -1, 0);

        final FermiLiteAssembly.Contig cleanContig = new FermiLiteAssembly.Contig(dummyContigSequence, dummyContigSequenceQuals, 100);
        final BwaMemAlignment cleanAlignment = new BwaMemAlignment(0, dummyRefId, 1000000, 1001000, 0, 1000, 60, 0, 100, 0, "1000M", "", "", -1, -1, 0);

        final FermiLiteAssembly.Contig contigWithGapInAlignment = new FermiLiteAssembly.Contig(dummyContigSequence, dummyContigSequenceQuals, 100);
        final BwaMemAlignment gappedAlignment = new BwaMemAlignment(0, dummyRefId,1000000, 1001000, 0, 1000, 60, 0, 100, 0, "700M50I200M50D50M", "", "", -1, -1, 0);

        final List<List<BwaMemAlignment>> allAlignments = Arrays.asList(Collections.singletonList(unmappedContigAlignment), Arrays.asList(firstAmbiguousMapping, secondAmbiguousMapping), Collections.singletonList(cleanAlignment), Collections.singletonList(gappedAlignment));
        final FermiLiteAssembly assembly = new FermiLiteAssembly(Arrays.asList(unmappedContig, contigWithAmbiguousMapping, cleanContig, contigWithGapInAlignment));

        final AlignedAssemblyOrExcuse collection = new AlignedAssemblyOrExcuse(1, assembly, allAlignments);

        final Iterable<Tuple2<Iterable<AlignmentRegion>, byte[]>> result = AssemblyAlignmentParser.forEachAssemblyNotExcuse(collection, refNames);
        final Iterator<Tuple2<Iterable<AlignmentRegion>, byte[]>> it = result.iterator();

        Assert.assertTrue(it.hasNext());
        Assert.assertFalse(it.next()._1().iterator().hasNext());
        Assert.assertFalse(it.next()._1().iterator().hasNext());

        final Iterator<AlignmentRegion> arOfCleanContigIt = it.next()._1().iterator();
        Assert.assertEquals(arOfCleanContigIt.next(),
                new AlignmentRegion("asm000001", "tig00002", new SimpleInterval(dummyRefName, 1000001, 1001000), TextCigarCodec.decode("1000M"), true, 60, 0, 1, 1000));
        Assert.assertFalse(arOfCleanContigIt.hasNext());

        Assert.assertEquals(Iterators.size(it.next()._1().iterator()), 3);

        Assert.assertEquals(AssemblyAlignmentParser.formatToAlignmentRegions(Collections.singleton(collection), refNames).size(), 2);
    }

    @Test
    public void testGappedAlignmentBreaker_OneInsertion() {

        final Cigar cigar = TextCigarCodec.decode("56S27M15I32M21S");
        final AlignmentRegion alignmentRegion = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 158), cigar, true, 60, 0, 57, 130);

        final List<AlignmentRegion> generatedARList = Utils.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, 1)).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.size(), 2);
        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 126), TextCigarCodec.decode("56S27M68S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM, 57, 83));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 127, 158), TextCigarCodec.decode("98S32M21S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM, 99, 130));
    }

    @Test
    public void testGappedAlignmentBreaker_OneDeletion() {
        final Cigar cigar = TextCigarCodec.decode("2S205M2D269M77S");
        final AlignmentRegion alignmentRegion = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 575), cigar, true, 60, 0, 208, 476);

        final List<AlignmentRegion> generatedARList = Utils.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, 1)).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.size(), 2);
        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 304), TextCigarCodec.decode("2S205M346S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM, 3, 207));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 307, 575), TextCigarCodec.decode("207S269M77S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM, 208, 476));
    }

    @Test
    public void testGappedAlignmentBreaker_Complex() {

        final AlignmentRegion alignmentRegion = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 414), TextCigarCodec.decode("397S118M2D26M6I50M7I26M1I8M13D72M398S"), true, 60, 65, 398, 711);

        final List<AlignmentRegion> generatedARList = Utils.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, 1)).collect(Collectors.toList());

        Assert.assertEquals(generatedARList.size(), 6);

        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 217), TextCigarCodec.decode("397S118M594S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM, 398, 515));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 220, 245), TextCigarCodec.decode("515S26M568S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM, 516, 541));
        Assert.assertEquals(generatedARList.get(2), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 246, 295), TextCigarCodec.decode("547S50M512S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM, 548, 597));
        Assert.assertEquals(generatedARList.get(3), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 296, 321), TextCigarCodec.decode("604S26M479S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM, 605, 630));
        Assert.assertEquals(generatedARList.get(4), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 322, 329), TextCigarCodec.decode("631S8M470S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM, 632, 639));
        Assert.assertEquals(generatedARList.get(5), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 343, 414), TextCigarCodec.decode("639S72M398S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM, 640, 711));
    }

    @Test
    public void testGappedAlignmentBreaker_GapSizeSensitivity() {

        final Cigar cigar = TextCigarCodec.decode("10M10D10M60I10M10I10M50D10M");
        final AlignmentRegion alignmentRegion = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 209), cigar, true, 60, 0, 1, 120);

        final List<AlignmentRegion> generatedARList = Utils.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, SVConstants.DiscoveryStepConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY)).collect(Collectors.toList());

        Assert.assertEquals(generatedARList.size(), 3);
        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 129), TextCigarCodec.decode("10M10D10M100S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM, 1, 20));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 130, 149), TextCigarCodec.decode("80S10M10I10M10S"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM, 81, 110));
        Assert.assertEquals(generatedARList.get(2), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 200, 209), TextCigarCodec.decode("110S10M"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM, 111, 120));
    }

    @Test
    public void testGappedAlignmentBreaker_HardAndSoftClip() {

        final AlignmentRegion alignmentRegion = new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 138), TextCigarCodec.decode("1H2S3M5I10M20D6M7S8H"), true, 60, 0, 4, 27);

        final List<AlignmentRegion> generatedARList = Utils.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, 1)).collect(Collectors.toList());

        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 100, 102), TextCigarCodec.decode("1H2S3M28S8H"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM, 4, 6));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 103, 112), TextCigarCodec.decode("1H10S10M13S8H"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM, 12, 21));
        Assert.assertEquals(generatedARList.get(2), new AlignmentRegion("1", "contig-1", new SimpleInterval("1", 133, 138), TextCigarCodec.decode("1H20S6M7S8H"), true, 60, SVConstants.DiscoveryStepConstants.MISSING_NM, 22, 27));
    }

    @Test
    public void testGappedAlignmentBreaker_TerminalInsertionOperatorToSoftClip() {

        // beginning with 'I'
        Cigar cigar = TextCigarCodec.decode("10I10M5I10M");
        AlignmentRegion alignmentRegion = new AlignmentRegion("1", "contig-1",
                new SimpleInterval("1", 101, 120), cigar,
                true, 60, 0, 11, 35);

        List<AlignmentRegion> generatedARList = Utils.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, 1)).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.size(), 2);
        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("1", "contig-1",
                new SimpleInterval("1", 101, 110), TextCigarCodec.decode("10S10M15S"),
                true, 60,
                SVConstants.DiscoveryStepConstants.MISSING_NM, 11, 20));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("1", "contig-1",
                new SimpleInterval("1", 111, 120), TextCigarCodec.decode("25S10M"),
                true, 60,
                SVConstants.DiscoveryStepConstants.MISSING_NM, 26, 35));

        // ending with 'I'
        cigar = TextCigarCodec.decode("10M5I10M10I");
        alignmentRegion = new AlignmentRegion("1", "contig-1",
                new SimpleInterval("1", 101, 120), cigar,
                true, 60, 0, 1, 25);

        generatedARList = Utils.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, 1)).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.size(), 2);
        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("1", "contig-1",
                new SimpleInterval("1", 101, 110), TextCigarCodec.decode("10M25S"),
                true, 60,
                SVConstants.DiscoveryStepConstants.MISSING_NM, 1, 10));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("1", "contig-1",
                new SimpleInterval("1", 111, 120), TextCigarCodec.decode("15S10M10S"),
                true, 60,
                SVConstants.DiscoveryStepConstants.MISSING_NM, 16, 25));
    }

    @Test
    public void testGappedAlignmentBreaker_TerminalInsertionNeighboringClippings(){

        Cigar cigar = TextCigarCodec.decode("10H20S30I40M50I60S70H");
        AlignmentRegion alignmentRegion = new AlignmentRegion("1", "contig-1",
                new SimpleInterval("1", 101, 140), cigar,
                true, 60, 0, 61, 100);
        List<AlignmentRegion> generatedARList = Utils.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, 1)).collect(Collectors.toList());
        // no internal gap, so nothing should change
        Assert.assertEquals(generatedARList.size(), 1);
        Assert.assertEquals(generatedARList.get(0), alignmentRegion);

        cigar = TextCigarCodec.decode("10H20S30I40M5D15M50I60S70H");
        alignmentRegion = new AlignmentRegion("1", "contig-1",
                new SimpleInterval("1", 101, 160), cigar,
                true, 60, 0, 61, 115);
        generatedARList = Utils.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, 1)).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.size(), 2);

        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("1", "contig-1",
                new SimpleInterval("1", 101, 140), TextCigarCodec.decode("10H50S40M125S70H"),
                true, 60,
                SVConstants.DiscoveryStepConstants.MISSING_NM, 61, 100));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("1", "contig-1",
                new SimpleInterval("1", 146, 160), TextCigarCodec.decode("10H90S15M110S70H"),
                true, 60,
                SVConstants.DiscoveryStepConstants.MISSING_NM, 101, 115));
    }

    @Test
    public void testGappedAlignmentBreaker_NegativeStrand() {
        // read data with AlignmentRegion.toString():
        // 19149	contig-8	chrUn_JTFH01000557v1_decoy	21	1459	-	10S1044M122I395M75I	60	11	1646	200
        final Cigar cigar = TextCigarCodec.decode("10S1044M122I395M75I");
        final AlignmentRegion alignmentRegion = new AlignmentRegion("19149", "contig-8",
                new SimpleInterval("chrUn_JTFH01000557v1_decoy", 21, 1459), cigar,
                false, 60, 200, 11, 1646);
        final List<AlignmentRegion> generatedARList = Utils.stream(AssemblyAlignmentParser.breakGappedAlignment(alignmentRegion, 1)).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.get(0), new AlignmentRegion("19149", "contig-8",
                new SimpleInterval("chrUn_JTFH01000557v1_decoy", 416, 1459), TextCigarCodec.decode("10S1044M592S"),
                false, 60,
                SVConstants.DiscoveryStepConstants.MISSING_NM, 11, 1054));
        Assert.assertEquals(generatedARList.get(1), new AlignmentRegion("19149", "contig-8",
                new SimpleInterval("chrUn_JTFH01000557v1_decoy", 21, 415), TextCigarCodec.decode("1176S395M75S"),
                false, 60,
                SVConstants.DiscoveryStepConstants.MISSING_NM, 1177, 1571));
    }

    @Test
    public void testGappedAlignmentBreaker_ExpectException() {
        int exceptionCount = 0;

        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10H10D10M"));} catch (final Exception ex) {++exceptionCount;}
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10S10D10M"));} catch (final Exception ex) {++exceptionCount;}
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10M10D10S"));} catch (final Exception ex) {++exceptionCount;}
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10M10D10H"));} catch (final Exception ex) {++exceptionCount;}

        try{willThrowOnInvalidCigar(TextCigarCodec.decode("10H10I10S"));} catch (final Exception ex) {++exceptionCount;}
        try{willThrowOnInvalidCigar(TextCigarCodec.decode("10H10D10S"));} catch (final Exception ex) {++exceptionCount;}

        try{willThrowOnInvalidCigar(TextCigarCodec.decode("10H10I10M10D10S"));} catch (final Exception ex) {++exceptionCount;}
        try{willThrowOnInvalidCigar(TextCigarCodec.decode("10H10D10M10I10S"));} catch (final Exception ex) {++exceptionCount;}

        // these 4 are fine now
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10H10I10M"));} catch (final Exception ex) {++exceptionCount;}
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10S10I10M"));} catch (final Exception ex) {++exceptionCount;}
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10M10I10S"));} catch (final Exception ex) {++exceptionCount;}
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10M10I10H"));} catch (final Exception ex) {++exceptionCount;}

        // last two are valid
        try{willThrowOnInvalidCigar(TextCigarCodec.decode("10H10M10I10M10S"));} catch (final Exception ex) {++exceptionCount;}
        try{willThrowOnInvalidCigar(TextCigarCodec.decode("10H10M10D10M10S"));} catch (final Exception ex) {++exceptionCount;}

        Assert.assertEquals(exceptionCount, 8);
    }

    private static Iterable<AlignmentRegion> willThrowOnInvalidCigar(final Cigar cigar) throws GATKException {
        final AlignmentRegion detailsDoesnotMatter = new AlignmentRegion("1", "contig-1",
                new SimpleInterval("1", 1, 110),
                cigar,
                true, 60, 0, 21, 30);
        return AssemblyAlignmentParser.breakGappedAlignment(detailsDoesnotMatter, 1);
    }

    @Test
    public void testCompactifyNeighboringSoftClippings() {
        Assert.assertEquals(new Cigar(AssemblyAlignmentParser.compactifyNeighboringSoftClippings(TextCigarCodec.decode("1H2S3S4M5D6M7I8M9S10S11H").getCigarElements())),
                TextCigarCodec.decode("1H5S4M5D6M7I8M19S11H"));
    }
}
