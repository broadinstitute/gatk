package org.broadinstitute.hellbender.tools.spark.sv.sga;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedAssembly;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SVDiscoveryTestDataProvider;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.AlignedAssemblyOrExcuse;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.MiniClusterUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.util.*;

public class SGATextFormatAlignmentParserUnitTest extends BaseTest {
    private static final String dummyRefName = "1";
    private static final int dummyRefId = Integer.valueOf(dummyRefName) - 1;
    private static final List<String> refNames = Collections.singletonList(dummyRefName);

    @Test(groups = "sv")
    @SuppressWarnings("deprecation")
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
            final Map<String, byte[]> contigNameAndSequence
                    = org.broadinstitute.hellbender.tools.spark.sv.sga.DiscoverVariantsFromContigAlignmentsSGASpark
                    .SGATextFormatAlignmentParser.extractContigNameAndSequenceFromTextFile(ctx, tempPath.toString()).collectAsMap();
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

            final List<AlignmentInterval> alignmentIntervalsForSimpleInversion = new ArrayList<>(8);
            final SimpleInterval referenceIntervalLeft = new SimpleInterval(refNames.get(0), alignmentStartsOnRef_0Based[2*pair]+1, alignmentStartsOnRef_0Based[2*pair]+cigars[2*pair].getReferenceLength()+1);
            final AlignmentInterval alignmentIntervalLeft = new AlignmentInterval(referenceIntervalLeft, alignmentStartsOnTig_0BasedInclusive[2*pair]+1, alignmentEndsOnTig_0BasedExclusive[2*pair],
                    strandedness[2*pair] ? cigars[2*pair] : CigarUtils.invertCigar(cigars[2*pair]),
                    strandedness[2*pair], mapQual[2*pair], mismatches[2*pair], 100, false, false);
            alignmentIntervalsForSimpleInversion.add(alignmentIntervalLeft);
            final SimpleInterval referenceIntervalRight = new SimpleInterval(refNames.get(0), alignmentStartsOnRef_0Based[2*pair+1]+1, alignmentStartsOnRef_0Based[2*pair+1]+cigars[2*pair+1].getReferenceLength()+1);
            final AlignmentInterval alignmentIntervalRight = new AlignmentInterval(referenceIntervalRight, alignmentStartsOnTig_0BasedInclusive[2*pair+1]+1, alignmentEndsOnTig_0BasedExclusive[2*pair+1],
                    strandedness[2*pair+1] ? cigars[2*pair+1] : CigarUtils.invertCigar(cigars[2*pair+1]),
                    strandedness[2*pair+1], mapQual[2*pair+1], mismatches[2*pair+1], 100, false, false);
            alignmentIntervalsForSimpleInversion.add(alignmentIntervalRight);

            if (pair == 0) {
                allContigs.add( new AlignedContig(AlignedAssemblyOrExcuse.formatContigName(0, 0), dummySequenceForContigOne, alignmentIntervalsForSimpleInversion, false) );
            } else if (pair <3) {
                allContigs.add( new AlignedContig(AlignedAssemblyOrExcuse.formatContigName(1, pair-1), pair==1 ? dummySequenceForContigTwo : dummySequenceForContigThree, alignmentIntervalsForSimpleInversion, false) );
            } else {
                allContigs.add( new AlignedContig(AlignedAssemblyOrExcuse.formatContigName(2, 0), dummySequenceForContigFour, alignmentIntervalsForSimpleInversion, false) );
            }
        }

        final Object[][] data = new Object[4][];
        data[0] = new Object[]{0, new AlignedAssembly(0, allContigs.subList(0, 1)), allContigs.subList(0, 1)};
        data[1] = new Object[]{1, new AlignedAssembly(1, allContigs.subList(1, 3)), allContigs.subList(1, 3)};
        data[2] = new Object[]{2, new AlignedAssembly(2, allContigs.subList(3, 4)), allContigs.subList(3, 4)};

        final List<AlignedContig> unmappedContig = Arrays.asList(new AlignedContig("asm000004:tig00001", SVDiscoveryTestDataProvider.makeDummySequence(20, (byte)'N'), Collections.emptyList(), false));
        data[3] = new Object[]{3, new AlignedAssembly(3, unmappedContig), unmappedContig};

        return data;
    }

    @Test(dataProvider = "AlignedAssemblyTextParserText", groups = "sv")
    @SuppressWarnings("deprecation")
    public void testEncodeAndDecodeAlignedAssemblyAsHadoopTextFileStringList(final Integer assemblyId, final AlignedAssembly expectedAssembly,
                                                                             final List<AlignedContig> contigs) {
        final Iterator<String> alignedContigStringIt
                = org.broadinstitute.hellbender.tools.spark.sv.sga.AlignAssembledContigsSpark.formatAlignedAssemblyAsText(expectedAssembly);
        if (assemblyId==1) {
            Assert.assertTrue(alignedContigStringIt.hasNext());
            final Tuple2<String, List<AlignmentInterval>> contigNameAndAlignments_0 =
                    org.broadinstitute.hellbender.tools.spark.sv.sga.DiscoverVariantsFromContigAlignmentsSGASpark
                            .SGATextFormatAlignmentParser.parseTextFileAlignmentIntervalLines(alignedContigStringIt.next());
            Assert.assertEquals(contigs.get(0).contigName, contigNameAndAlignments_0._1());
            Assert.assertEquals(contigs.get(0).alignmentIntervals, contigNameAndAlignments_0._2());
            final Tuple2<String, List<AlignmentInterval>> contigNameAndAlignments_1 =
                    org.broadinstitute.hellbender.tools.spark.sv.sga.DiscoverVariantsFromContigAlignmentsSGASpark
                            .SGATextFormatAlignmentParser.parseTextFileAlignmentIntervalLines(alignedContigStringIt.next());
            Assert.assertEquals(contigs.get(1).contigName, contigNameAndAlignments_1._1());
            Assert.assertEquals(contigs.get(1).alignmentIntervals, contigNameAndAlignments_1._2());
        } else {
            Assert.assertTrue(alignedContigStringIt.hasNext());
            final Tuple2<String, List<AlignmentInterval>> contigNameAndAlignments =
                    org.broadinstitute.hellbender.tools.spark.sv.sga.DiscoverVariantsFromContigAlignmentsSGASpark
                            .SGATextFormatAlignmentParser.parseTextFileAlignmentIntervalLines(alignedContigStringIt.next());
            Assert.assertEquals(contigs.get(0).contigName, contigNameAndAlignments._1());
            Assert.assertEquals(contigs.get(0).alignmentIntervals, contigNameAndAlignments._2());
        }
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
    @SuppressWarnings("deprecation")
    public void testEnCodeAndDecodeSimpleIntervalAsString_valid(final SimpleInterval interval) {
        Assert.assertEquals(
                org.broadinstitute.hellbender.tools.spark.sv.sga.AlignAssembledContigsSpark.decodeStringAsSimpleInterval(
                        org.broadinstitute.hellbender.tools.spark.sv.sga.AlignAssembledContigsSpark.encodeSimpleIntervalAsString(interval)),
                interval);
    }
}
