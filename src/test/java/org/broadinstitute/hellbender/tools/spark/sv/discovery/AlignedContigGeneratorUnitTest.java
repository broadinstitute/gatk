package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.collect.Iterables;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.SVConstants;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryPipelineSpark;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.AlignedAssemblyOrExcuse;
import org.broadinstitute.hellbender.tools.spark.sv.sga.AlignAssembledContigsSpark;
import org.broadinstitute.hellbender.tools.spark.sv.sga.DiscoverVariantsFromContigAlignmentsSGASpark;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.MiniClusterUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.testng.Assert.assertEquals;

public class AlignedContigGeneratorUnitTest extends BaseTest {
    private static final String dummyRefName = "1";
    private static final int dummyRefId = Integer.valueOf(dummyRefName) - 1;
    private static final List<String> refNames = Collections.singletonList(dummyRefName);

    @Test(groups = "sv")
    public void testExtractContigNameAndSequenceFromTextFile() throws Exception {
        MiniClusterUtils.runOnIsolatedMiniCluster(cluster -> {
            //use the HDFS on the mini cluster
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);
            final Path tempPath = new Path(workingDirectory, "testLocalAssemblyFormatRead");
            final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

            final FileSystem fs = tempPath.getFileSystem(ctx.hadoopConfiguration());
            final FSDataOutputStream fsOutStream = fs.create(tempPath);

            final String packedFasta = ">contig-4 521 0|ATTTGTGAGCGCTTTGAGGCCTTTGGAGAAAAAGGAAATGTATTCACAGAAAAACTTGAAAAAAAGCTTCTGGTAAACTGT" +
                    "TTTGTAATGTGTACAATCATTTCACAGAGATCAGTGTTTCTTTTCATTGAGCAGCTTGGAAACTCTATTGTTGTAGAATCTGCAAACGGATATTTTTCAGTG" +
                    "CTTTGATGCCTGTGTTAAAAAAGGAAATATCCTCACATAAAAAATGACAGAAAATTTATGAGAAACTTCTTTGTGATGTGTGCATTTATGCCACAGAATTGA" +
                    "ACCATTTTATGATTGAGCAGTTTGGAAACAGTCTTTTTGTGGAATCTAAAAAGAGATATTTATGAGCGCATTGAGGCCTACAGTAAAAAAGGAAATATCTTC" +
                    "ACATAAAAACTAGCAAGGAGCTTTCTAAGAAACTGCTTTGTGATGCGTGAATTCATCTCACAGAGGTAAATGTTTCTTTGCATTGAACAGTGGAAACTCTGT" +
                    "TCTTGTAGAATCTGCAAAGTGATATTTGTGAG|>contig-11 164 0|TCACAGAATTGAACGACCTCTTTGTTTGAGCAGTTTGGGAACAGCCTTTTTG" +
                    "TAGATTCTGCAAAGGGATATTTGTAAGCCCTTTGAGGACTATGGTGAAAACGTAAATATCTTCACATAACTAGACAGAAGGTTGCTGAAAAGCTGCTTTGTG" +
                    "ATGTGTGATT|>contig-0 207 0|GCATTGAACAGTGAAAACTCTGTTCTTGTAGAATCTGCAAAGTGATATTTGTGAGTGTTTTGAGGCCTATGGTGA" +
                    "AAAAGGAAATATCTTCAGAAAAACTAGACAGAAGCTTTCTGAGAATATTCTTTGTGATATATGCATTCATCTCACAGAATTGAACGACCTCTTTGTTTGAGC" +
                    "AGTTTGGGAACAGCCTTTTTGTAGATTCTG";
            final String oneAssembly = "1\t" + packedFasta;

            fsOutStream.writeBytes(oneAssembly);
            fsOutStream.writeBytes("\n");
            fsOutStream.close();
            fs.deleteOnExit(tempPath);

            Assert.assertTrue(SparkUtils.pathExists(ctx, tempPath));
            final Map<String, byte[]> contigNameAndSequence = DiscoverVariantsFromContigAlignmentsSGASpark.SGATextFormatAlignmentParser.extractContigNameAndSequenceFromTextFile(ctx, tempPath.toString()).collectAsMap();
            Assert.assertEquals(contigNameAndSequence.size(), 3);
            Assert.assertEquals(contigNameAndSequence.get("asm000001:tig00000"), "GCATTGAACAGTGAAAACTCTGTTCTTGTAGAATCTGCAAAGTGATATTTGTGAGTGTTTTGAGGCCTATGGTGAAAAAGGAAATATCTTCAGAAAAACTAGACAGAAGCTTTCTGAGAATATTCTTTGTGATATATGCATTCATCTCACAGAATTGAACGACCTCTTTGTTTGAGCAGTTTGGGAACAGCCTTTTTGTAGATTCTG".getBytes());
            Assert.assertEquals(contigNameAndSequence.get("asm000001:tig00004"), "ATTTGTGAGCGCTTTGAGGCCTTTGGAGAAAAAGGAAATGTATTCACAGAAAAACTTGAAAAAAAGCTTCTGGTAAACTGTTTTGTAATGTGTACAATCATTTCACAGAGATCAGTGTTTCTTTTCATTGAGCAGCTTGGAAACTCTATTGTTGTAGAATCTGCAAACGGATATTTTTCAGTGCTTTGATGCCTGTGTTAAAAAAGGAAATATCCTCACATAAAAAATGACAGAAAATTTATGAGAAACTTCTTTGTGATGTGTGCATTTATGCCACAGAATTGAACCATTTTATGATTGAGCAGTTTGGAAACAGTCTTTTTGTGGAATCTAAAAAGAGATATTTATGAGCGCATTGAGGCCTACAGTAAAAAAGGAAATATCTTCACATAAAAACTAGCAAGGAGCTTTCTAAGAAACTGCTTTGTGATGCGTGAATTCATCTCACAGAGGTAAATGTTTCTTTGCATTGAACAGTGGAAACTCTGTTCTTGTAGAATCTGCAAAGTGATATTTGTGAG".getBytes());
            Assert.assertEquals(contigNameAndSequence.get("asm000001:tig00011"), "TCACAGAATTGAACGACCTCTTTGTTTGAGCAGTTTGGGAACAGCCTTTTTGTAGATTCTGCAAAGGGATATTTGTAAGCCCTTTGAGGACTATGGTGAAAACGTAAATATCTTCACATAACTAGACAGAAGGTTGCTGAAAAGCTGCTTTGTGATGTGTGATT".getBytes());
        });
    }

    /**
     * 4 contigs of 3 assemblies pointing to the same inversion event.
     *  contig 1 belongs to assembly 1
     *  contig 2 and 3 belong to assembly 2
     *  contig 4 belongs to assembly 3.
     */
    @DataProvider(name = "AlignedAssemblyTextParserText")
    private Object[][] createInputsAndExpectedResults() {

        final int[] alignmentStartsOnRef_0Based = {96, 196, 195, 95, 101, 201, 101, 201};
        final int[] alignmentStartsOnTig_0BasedInclusive = {0, 4, 0, 5, 0, 6, 0, 7};
        final int[] alignmentEndsOnTig_0BasedExclusive = {4, 8, 5, 10, 6, 12, 7, 14};
        final int[] seqLen = {8, 8, 10, 10, 12, 12, 14, 14};
        final int[] mapQual = {0, 1, 10, 20, 30, 40, 50, 60};
        final int[] mismatches = {0, 1, 1, 0, 2, 3, 3, 2};
        final boolean[] strandedness = {true, false, true, false, false, true, false, true};
        final String[] cigarStrings = {"4M4S", "4M4H", "5M5S", "5M5H", "6S6M", "6H6M", "7S7M", "7H7M"}; // each different number represent a different contig's pair of chimeric alignments
        final Cigar[] cigars = Arrays.stream(cigarStrings).map(TextCigarCodec::decode).toArray(Cigar[]::new);

        // these sequence are technically wrong the for the inversion event, but the test purpose is for serialization so it is irrelevant
        final byte[] dummySequenceForContigOne = SVDiscoveryTestDataProvider.makeDummySequence(seqLen[0], (byte)'A');

        final byte[] dummySequenceForContigTwo = SVDiscoveryTestDataProvider.makeDummySequence(seqLen[0], (byte)'T');
        final byte[] dummySequenceForContigThree = SVDiscoveryTestDataProvider.makeDummySequence(seqLen[0], (byte)'C');

        final byte[] dummySequenceForContigFour = SVDiscoveryTestDataProvider.makeDummySequence(seqLen[0], (byte)'G');


        final List<AlignedContig> allContigs = new ArrayList<>();

        for(int pair=0; pair<cigars.length/2; ++pair) {

            final List<AlignedAssembly.AlignmentInterval> alignmentIntervalsForSimpleInversion = new ArrayList<>(8);
            final SimpleInterval referenceIntervalLeft = new SimpleInterval(refNames.get(0), alignmentStartsOnRef_0Based[2*pair]+1, alignmentStartsOnRef_0Based[2*pair]+cigars[2*pair].getReferenceLength()+1);
            final AlignedAssembly.AlignmentInterval alignmentIntervalLeft = new AlignedAssembly.AlignmentInterval(referenceIntervalLeft, alignmentStartsOnTig_0BasedInclusive[2*pair]+1, alignmentEndsOnTig_0BasedExclusive[2*pair],
                    strandedness[2*pair] ? cigars[2*pair] : CigarUtils.invertCigar(cigars[2*pair]),
                    strandedness[2*pair], mapQual[2*pair], mismatches[2*pair]);
            alignmentIntervalsForSimpleInversion.add(alignmentIntervalLeft);
            final SimpleInterval referenceIntervalRight = new SimpleInterval(refNames.get(0), alignmentStartsOnRef_0Based[2*pair+1]+1, alignmentStartsOnRef_0Based[2*pair+1]+cigars[2*pair+1].getReferenceLength()+1);
            final AlignedAssembly.AlignmentInterval alignmentIntervalRight = new AlignedAssembly.AlignmentInterval(referenceIntervalRight, alignmentStartsOnTig_0BasedInclusive[2*pair+1]+1, alignmentEndsOnTig_0BasedExclusive[2*pair+1],
                    strandedness[2*pair+1] ? cigars[2*pair+1] : CigarUtils.invertCigar(cigars[2*pair+1]),
                    strandedness[2*pair+1], mapQual[2*pair+1], mismatches[2*pair+1]);
            alignmentIntervalsForSimpleInversion.add(alignmentIntervalRight);

            if (pair == 0) {
                allContigs.add( new AlignedContig(AlignedAssemblyOrExcuse.formatContigName(0, 0), dummySequenceForContigOne, alignmentIntervalsForSimpleInversion) );
            } else if (pair <3) {
                allContigs.add( new AlignedContig(AlignedAssemblyOrExcuse.formatContigName(1, pair-1), pair==1 ? dummySequenceForContigTwo : dummySequenceForContigThree, alignmentIntervalsForSimpleInversion) );
            } else {
                allContigs.add( new AlignedContig(AlignedAssemblyOrExcuse.formatContigName(2, 0), dummySequenceForContigFour, alignmentIntervalsForSimpleInversion) );
            }
        }

        final Object[][] data = new Object[4][];
        data[0] = new Object[]{0, new AlignedAssembly(0, allContigs.subList(0, 1)), allContigs.subList(0, 1)};
        data[1] = new Object[]{1, new AlignedAssembly(1, allContigs.subList(1, 3)), allContigs.subList(1, 3)};
        data[2] = new Object[]{2, new AlignedAssembly(2, allContigs.subList(3, 4)), allContigs.subList(3, 4)};

        final List<AlignedContig> unmappedContig = Arrays.asList(new AlignedContig("asm000004:tig00001", SVDiscoveryTestDataProvider.makeDummySequence(20, (byte)'N'), Collections.emptyList()));
        data[3] = new Object[]{3, new AlignedAssembly(3, unmappedContig), unmappedContig};

        return data;
    }

    @Test(dataProvider = "AlignedAssemblyTextParserText", groups = "sv")
    public void testEncodeAndDecodeAlignedAssemblyAsHadoopTextFileStringList(final Integer assemblyId, final AlignedAssembly expectedAssembly,
                                                                             final List<AlignedContig> contigs) {
        final Iterator<String> alignedContigStringIt = AlignAssembledContigsSpark.formatAlignedAssemblyAsText(expectedAssembly);
        if (assemblyId==1) {
            Assert.assertTrue(alignedContigStringIt.hasNext());
            final Tuple2<String, List<AlignedAssembly.AlignmentInterval>> contigNameAndAlignments_0 =
                    DiscoverVariantsFromContigAlignmentsSGASpark.SGATextFormatAlignmentParser.parseTextFileAlignmentIntervalLines(alignedContigStringIt.next());
            Assert.assertEquals(contigs.get(0).contigName, contigNameAndAlignments_0._1());
            Assert.assertEquals(contigs.get(0).alignmentIntervals, contigNameAndAlignments_0._2());
            final Tuple2<String, List<AlignedAssembly.AlignmentInterval>> contigNameAndAlignments_1 =
                    DiscoverVariantsFromContigAlignmentsSGASpark.SGATextFormatAlignmentParser.parseTextFileAlignmentIntervalLines(alignedContigStringIt.next());
            Assert.assertEquals(contigs.get(1).contigName, contigNameAndAlignments_1._1());
            Assert.assertEquals(contigs.get(1).alignmentIntervals, contigNameAndAlignments_1._2());
        } else {
            Assert.assertTrue(alignedContigStringIt.hasNext());
            final Tuple2<String, List<AlignedAssembly.AlignmentInterval>> contigNameAndAlignments =
                    DiscoverVariantsFromContigAlignmentsSGASpark.SGATextFormatAlignmentParser.parseTextFileAlignmentIntervalLines(alignedContigStringIt.next());
            Assert.assertEquals(contigs.get(0).contigName, contigNameAndAlignments._1());
            Assert.assertEquals(contigs.get(0).alignmentIntervals, contigNameAndAlignments._2());
        }
    }

    @Test(groups = "sv")
    public void testParseAlignedAssembledContigTextLine() throws Exception {
        MiniClusterUtils.runOnIsolatedMiniCluster(cluster -> {
            //use the HDFS on the mini cluster
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);
            final Path tempPath = new Path(workingDirectory, "testLocalAssemblyFormatRead");
            final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

            final FileSystem fs = tempPath.getFileSystem(ctx.hadoopConfiguration());
            final FSDataOutputStream fsOutStream = fs.create(tempPath);

            final String gappedAlignmentContig_1 = "asm000001:tig00001\t1-200%CTG=1START=101END=300%45M100D55M%+%60%3";
            final String gappedAlignmentContig_2 = "asm000001:tig00002\t1-200%CTG=1START=106END=305%60M100D40M%-%60%5";
            fsOutStream.writeBytes(gappedAlignmentContig_1);
            fsOutStream.writeBytes("\n");
            fsOutStream.writeBytes(gappedAlignmentContig_2);
            fsOutStream.writeBytes("\n");
            fsOutStream.close();
            fs.deleteOnExit(tempPath);

            Assert.assertTrue(SparkUtils.pathExists(ctx, tempPath));

            final Map<String, List<AlignedAssembly.AlignmentInterval>> contigNameAndAlignments
                    = DiscoverVariantsFromContigAlignmentsSGASpark.SGATextFormatAlignmentParser.parseAndBreakAlignmentTextRecords(ctx, tempPath.toString(), null).collectAsMap();

            Assert.assertEquals(contigNameAndAlignments.keySet().size(), 2);
            Assert.assertEquals(contigNameAndAlignments.get("asm000001:tig00001").size(), 2);
            Assert.assertEquals(contigNameAndAlignments.get("asm000001:tig00002").size(), 2);
        });
    }

    @Test(groups = "sv")
    public void testConvertGATKReadsToAlignedContig() throws Exception {
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

        final List<SAMRecord> reads = Stream.of(read1, read2, read3).map(read -> read.convertToSAMRecord(null)).collect(Collectors.toList());

        final AlignedContig alignedContig = DiscoverVariantsFromContigAlignmentsSAMSpark.SAMFormattedContigAlignmentParser.parseReadsAndBreakGaps(reads, null, SVConstants.DiscoveryStepConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, null);
        assertEquals(alignedContig.contigSequence, read2.getBases());

        assertEquals(alignedContig.alignmentIntervals.size(), 3);

        final byte[] read4Bytes = SVDiscoveryTestDataProvider.LONG_CONTIG1.getBytes();
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

        List<SAMRecord> reads2 = Stream.of(read4, read5).map(read -> read.convertToSAMRecord(null)).collect(Collectors.toList());

        final AlignedContig alignedContig2 = DiscoverVariantsFromContigAlignmentsSAMSpark.SAMFormattedContigAlignmentParser.parseReadsAndBreakGaps(reads2, null, SVConstants.DiscoveryStepConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, null);
        // these should be the reverse complements of each other
        assertEquals(alignedContig2.contigSequence.length, read4.getBases().length);

        final List<AlignedAssembly.AlignmentInterval> alignmentIntervals2 = alignedContig2.alignmentIntervals;
        assertEquals(alignmentIntervals2.size(), 2);

        final AlignedAssembly.AlignmentInterval alignmentInterval4 = alignmentIntervals2.get(0);
        assertEquals(alignmentInterval4.startInAssembledContig, 1);
        assertEquals(alignmentInterval4.endInAssembledContig, read4Seq.length() - 1986);

        final AlignedAssembly.AlignmentInterval alignmentInterval5 = alignmentIntervals2.get(1);
        assertEquals(alignmentInterval5.startInAssembledContig, 3604);
        assertEquals(alignmentInterval5.endInAssembledContig, read4Seq.length());
    }

    @Test(groups = "sv")
    public void testConvertAlignedAssemblyOrExcuseToAlignedContigsDirectAndConcordanceWithSAMRoute() {

        // test "failed" assembly doesn't produce anything
        final AlignedAssemblyOrExcuse excuse = new AlignedAssemblyOrExcuse(1, "justATest");
        Assert.assertTrue(StructuralVariationDiscoveryPipelineSpark.InMemoryAlignmentParser.filterAndConvertToAlignedContigDirect(Collections.singletonList(excuse), refNames, null).isEmpty());

        // produce test assembly and alignment
        final byte[] dummyContigSequence = SVDiscoveryTestDataProvider.makeDummySequence(1000, (byte)'T');
        final byte[] dummyContigSequenceQuals = SVDiscoveryTestDataProvider.makeDummySequence(1000, (byte)'A');
        final List<FermiLiteAssembly.Connection> dummyConnections = Collections.emptyList();

        final FermiLiteAssembly.Contig unmappedContig = new FermiLiteAssembly.Contig(dummyContigSequence, dummyContigSequenceQuals, 100); // totally random 100 supporting reads
        unmappedContig.setConnections(dummyConnections);
        final BwaMemAlignment unmappedContigAlignment = new BwaMemAlignment(4, -1, -1, -1, -1, -1, -1, -1, 0, 0, "", "", "", -1, -1, 0);

        final FermiLiteAssembly.Contig contigWithAmbiguousMapping = new FermiLiteAssembly.Contig(dummyContigSequence, dummyContigSequenceQuals, 100);
        contigWithAmbiguousMapping.setConnections(dummyConnections);
        final BwaMemAlignment firstAmbiguousMapping = new BwaMemAlignment(256, dummyRefId, 1000000, 1001000, 0, 1000, 0, 20, 100, 100, "800M50I100M50D50M", "", "", -1, -1, 0); // technically not correct but doesn't matter for this case
        final BwaMemAlignment secondAmbiguousMapping = new BwaMemAlignment(272, dummyRefId, 2000000, 2001000, 0, 1000, 0, 50, 100, 100, "700M50I200M50D50M", "", "", -1, -1, 0);

        final FermiLiteAssembly.Contig cleanContig = new FermiLiteAssembly.Contig(dummyContigSequence, dummyContigSequenceQuals, 100);
        cleanContig.setConnections(dummyConnections);
        final BwaMemAlignment cleanAlignment = new BwaMemAlignment(0, dummyRefId, 1000000, 1001000, 0, 1000, 60, 0, 100, 0, "1000M", "", "", -1, -1, 0);

        final FermiLiteAssembly.Contig contigWithGapInAlignment = new FermiLiteAssembly.Contig(dummyContigSequence, dummyContigSequenceQuals, 100);
        contigWithGapInAlignment.setConnections(dummyConnections);
        final BwaMemAlignment gappedAlignment = new BwaMemAlignment(0, dummyRefId,1000000, 1001000, 0, 1000, 60, 0, 100, 0, "700M50I200M50D50M", "", "", -1, -1, 0);

        final List<List<BwaMemAlignment>> allAlignments = Arrays.asList(Collections.singletonList(unmappedContigAlignment), Arrays.asList(firstAmbiguousMapping, secondAmbiguousMapping), Collections.singletonList(cleanAlignment), Collections.singletonList(gappedAlignment));
        final FermiLiteAssembly assembly = new FermiLiteAssembly(Arrays.asList(unmappedContig, contigWithAmbiguousMapping, cleanContig, contigWithGapInAlignment));
        final AlignedAssemblyOrExcuse alignedAssembly = new AlignedAssemblyOrExcuse(1, assembly, allAlignments);

        // test contig extraction without unmapped and unambiguous filtering
        final Iterable<AlignedContig> alignedContigsIncludingUnmapped = StructuralVariationDiscoveryPipelineSpark.InMemoryAlignmentParser.getAlignedContigsInOneAssembly(alignedAssembly, refNames, null);
        Assert.assertEquals(Iterables.size(alignedContigsIncludingUnmapped), 4);

        final Iterator<AlignedContig> it = alignedContigsIncludingUnmapped.iterator();

        Assert.assertTrue(it.next().alignmentIntervals.isEmpty());
        Assert.assertTrue(it.next().alignmentIntervals.isEmpty());

        final List<AlignedAssembly.AlignmentInterval> alignmentIntervalsForCleanContig = it.next().alignmentIntervals;
        Assert.assertEquals(alignmentIntervalsForCleanContig.size(), 1);
        Assert.assertEquals(alignmentIntervalsForCleanContig.get(0), new AlignedAssembly.AlignmentInterval(new SimpleInterval(dummyRefName, 1000001, 1001000), 1, 1000, TextCigarCodec.decode("1000M"), true, 60, 0));

        final List<AlignedAssembly.AlignmentInterval> alignmentIntervalsForContigWithGappedAlignment = it.next().alignmentIntervals;
        Assert.assertEquals(alignmentIntervalsForContigWithGappedAlignment.size(), 3);

        // test direct conversion (essentially the filtering step)
        final List<AlignedContig> parsedContigsViaDirectRoute
                = StructuralVariationDiscoveryPipelineSpark.InMemoryAlignmentParser.filterAndConvertToAlignedContigDirect(Collections.singleton(alignedAssembly), refNames, null);
        Assert.assertEquals(parsedContigsViaDirectRoute.size(), 2);
        Assert.assertTrue( parsedContigsViaDirectRoute.containsAll(Utils.stream(alignedContigsIncludingUnmapped).filter(ctg -> !ctg.alignmentIntervals.isEmpty()).collect(Collectors.toList())) );

        // concordance test with results obtained via SAM route
        final List<AlignedContig> parsedContigsViaSAMRoute
                = StructuralVariationDiscoveryPipelineSpark.InMemoryAlignmentParser.filterAndConvertToAlignedContigViaSAM(Collections.singletonList(alignedAssembly), hg19Header, SparkContextFactory.getTestSparkContext(), null).collect();
        Assert.assertEquals(parsedContigsViaDirectRoute, parsedContigsViaSAMRoute);
    }

    @DataProvider(name = "InvalidSimpleIntervalStrings")
    private Object[][] createInvalidSimpleIntervalStrings() {
        return new Object[][]{
                {"fjdskjfklsdjf"},
                {"1START=20000END=20010"},
                {"CTG=1START=20000END=2d0"},
                {"CTG=1START=20fd0END=200000"},
                {"CTG=1START=20END=19"}
        };
    }

    @Test(dataProvider = "InvalidSimpleIntervalStrings", expectedExceptions = GATKException.class)
    public void testEnCodeAndDecodeSimpleIntervalAsString_expectException(final String text) {
        AlignAssembledContigsSpark.decodeStringAsSimpleInterval(text);
    }

    @DataProvider(name = "ValidSimpleIntervalStrings")
    private Object[][] createSimpleIntervalStrings() {
        final SimpleInterval first = new SimpleInterval("1", 20000, 20010);
        final SimpleInterval second = new SimpleInterval("chr1", 20000, 20010);
        final SimpleInterval third = new SimpleInterval("HLA-DRB1*15:03:01:02", 20000, 20010);
        final SimpleInterval fourth = new SimpleInterval("chrUn_JTFH01001998v1_decoy", 20000, 20010);
        final SimpleInterval fifth = new SimpleInterval("chr19_KI270938v1_alt", 20000, 20010);
        final SimpleInterval sixth = new SimpleInterval("chrEBV", 20000, 20010);
        final SimpleInterval seventh = new SimpleInterval("chrUn_KI270529v1", 20000, 20010);
        final SimpleInterval eighth = new SimpleInterval("chr1_KI270713v1_random", 20000, 20010);
        final SimpleInterval ninth = new SimpleInterval("chrM", 20000, 20010);
        final SimpleInterval tenth = new SimpleInterval("chrX", 20000, 20010);

        return new Object[][]{{first}, {second} , {third}, {fourth}, {fifth}, {sixth}, {seventh}, {eighth}, {ninth}, {tenth}};
    }

    @Test(dataProvider = "ValidSimpleIntervalStrings")
    public void testEnCodeAndDecodeSimpleIntervalAsString_valid(final SimpleInterval interval) {
        Assert.assertEquals(AlignAssembledContigsSpark.decodeStringAsSimpleInterval(AlignAssembledContigsSpark.encodeSimpleIntervalAsString(interval)), interval);
    }
}
