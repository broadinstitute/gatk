package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqRecord;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.junit.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

/**
 * Created by shuang on 9/19/16.
 */
public final class SVFastqUtilsUnitTest extends BaseTest{

    private static JavaSparkContext ctx;

    private static GATKRead testMappedRead;
    private static FastqRecord testMappedReadFastq;

    private static GATKRead testUnmappedRead;
    private static FastqRecord testUnmappedReadFastq;

    private static GATKRead testMultiMappedNonprimaryRead;
    private static FastqRecord testMultiMappedNonprimaryReadFastq;

    @BeforeClass
    private static void setupSparkAndTestRead(){
        SparkContextFactory.enableTestSparkContext();
        ctx = SparkContextFactory.getTestSparkContext(Collections.emptyMap());
        generateTestReads();
    }

    static void generateTestReads(){
        //Using the following actual read
        /*"HJYFJCCXX160204:5:2214:27793:55719 83  20  799852  60  4S77M4I30M1D3M1I15M6D17M    =   799624  -377"
        "CTTGCGGTGGGAGAGGGGAACCTGCGCCCGCGATTCTCTCTCTCTCCCCCTCTCTCCCTCTCTCTTTCTCGCCCCCCCCCTCCCCCTCTCTCTCCCTCCCTCCCTCTCTCTTTCTCTCTCCCCCTCCCTCCCTCTCTCTCCCTCCCTCCCT"
        "######################@9-5;4%=0'..=).(<:=).7.?<>8>>).3,=3).=<).2..:).%$963:=='&9<)&93;.989.&.6='.9.:=>49,980<=>.<9.1;98====9,-17/==<9=5==>=>=>=94>?>9<>"
        "MC:Z:151M  MD:Z:8A15T41T7T27C1C2^A4T13^TCTCTG17    PG:Z:MarkDuplicates.1I  RG:Z:HJYFJ.5    NM:i:19 MQ:i:60"
        "OQ:Z:######################FA7<A7,F<,,,<,,,FAA,,<(KFKAKF,,7,F7,,FF,,7,,A,(,(AA7AFF,(AA,(A7F,A7A,,,AF,,A7FFKAA,A7<AFF,KA,A<A7KKFFA,<7A,KFF<F<FKKAKFKKA<FFFAAA"
        "UQ:i:148   AS:i:72"*/
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        byte[] seq = "CTTGCGGTGGGAGAGGGGAACCTGCGCCCGCGATTCTCTCTCTCTCCCCCTCTCTCCCTCTCTCTTTCTCGCCCCCCCCCTCCCCCTCTCTCTCCCTCCCTCCCTCTCTCTTTCTCTCTCCCCCTCCCTCCCTCTCTCTCCCTCCCTCCCT".getBytes();
        byte[] qual = SAMUtils.fastqToPhred("######################@9-5;4%=0'..=).(<:=).7.?<>8>>).3,=3).=<).2..:).%$963:=='&9<)&93;.989.&.6='.9.:=>49,980<=>.<9.1;98====9,-17/==<9=5==>=>=>=94>?>9<>");
        testMappedRead = ArtificialReadUtils.createArtificialRead(header, "HJYFJCCXX160204:5:2214:27793:55719", "20", 799852, seq, qual, "4S77M4I30M1D3M1I15M6D17M");
        testMappedRead.setIsFirstOfPair();
        testMappedReadFastq = new FastqRecord("@HJYFJCCXX160204:5:2214:27793:55719/1 mapping=20:799852;4S77M4I30M1D3M1I15M6D17M", new String(seq), "+", SAMUtils.phredToFastq(qual));

        /*"HJYFJCCXX160204:7:2103:31669:12754	133	20	700612	0	*	=	700612	0"
        "GGCATTTGTAGTGCCCTGACCGACGCCACCCACTGCTCTGCAGCCACATATCATGGCACTCCACAGGCAGCAGCCCAGGGAGCCAGTGCCAGCGTTCTCCGCGTTTCACCCCCGCCACAGGGCCGCCCCGCGCCGCCGCGTTCCCCGCCCC"
        "<*7>??>4>+)>===%=>7:###################################################################################################################################"
        "MC:Z:151M	PG:Z:MarkDuplicates.1K	RG:Z:HJYFJ.7"
        "OQ:Z:A,AAFKK<F,,FKKK(AF7F###################################################################################################################################"*/
        seq = "GGCATTTGTAGTGCCCTGACCGACGCCACCCACTGCTCTGCAGCCACATATCATGGCACTCCACAGGCAGCAGCCCAGGGAGCCAGTGCCAGCGTTCTCCGCGTTTCACCCCCGCCACAGGGCCGCCCCGCGCCGCCGCGTTCCCCGCCCC".getBytes();
        qual = SAMUtils.fastqToPhred("<*7>??>4>+)>===%=>7:###################################################################################################################################");
        testUnmappedRead = ArtificialReadUtils.createArtificialUnmappedReadWithAssignedPosition(header, "20", 700612, seq, qual);
        testUnmappedRead.setIsFirstOfPair();
        testUnmappedReadFastq = new FastqRecord("@HJYFJCCXX160204:7:2103:31669:12754/1", new String(seq), "+", new String(qual));

        /*"HJYFJCCXX160204:6:2219:2686:64404	329	20	799921	0	55S46M2I21M4I16M7S	=	799921	0"
        "CCTCTCTCTCTCTCTGCATCCCTCCCTCTCTCTCTCTCTCTCTGCCTCCCTCCCTCCCCCTCTCTCTCTCTCCCTCCCTCCCTCTCTCTCCCTCCCCCCCCCTCTCCCTCCCTCTCTCTCTCTCCCNTTCCCTCCCTCCCTCCCCCCCTCC"
        "<:?9@>?=?=?=>=?=<>=<=<<<<=><><><><?=?=>=>=?>>7?8;3<;3605$7%4=<:;0;<:<801(1?4&7<<%;;;6:==?##############################################################"
        "SA:Z:2,1723980,-,62S89M,49,4;	MD:Z:33C0T3A1T2T18G20	PG:Z:MarkDuplicates.1H	RG:Z:HJYFJ.6	NM:i:12"
        "OQ:Z:<AFAFKKKKKKKKKKKFKKFKFAKFKKKKKKKKKKKKKFKKKKKKAKAF<<F<A,<(A(<AF7F,F<F<A,<,7FA(A7F(F7F7FFKK##############################################################"
        "UQ:i:52	AS:i:35"*/
        seq = "CCTCTCTCTCTCTCTGCATCCCTCCCTCTCTCTCTCTCTCTCTGCCTCCCTCCCTCCCCCTCTCTCTCTCTCCCTCCCTCCCTCTCTCTCCCTCCCCCCCCCTCTCCCTCCCTCTCTCTCTCTCCCNTTCCCTCCCTCCCTCCCCCCCTCC".getBytes();
        qual = SAMUtils.fastqToPhred("<:?9@>?=?=?=>=?=<>=<=<<<<=><><><><?=?=>=>=?>>7?8;3<;3605$7%4=<:;0;<:<801(1?4&7<<%;;;6:==?##############################################################");
        testMultiMappedNonprimaryRead = ArtificialReadUtils.createArtificialRead(header, "HJYFJCCXX160204:6:2219:2686:64404", "20", 799921, seq, qual, "55S46M2I21M4I16M7S");
        testMultiMappedNonprimaryRead.setIsSupplementaryAlignment(true);
        testMultiMappedNonprimaryRead.setIsFirstOfPair();
        testMultiMappedNonprimaryReadFastq = new FastqRecord("@HJYFJCCXX160204:6:2219:2686:64404/1", new String(seq), "+", new String(qual));
    }

    @AfterClass
    private static void closeSpark(){
        SparkContextFactory.stopSparkContext(ctx);
    }

    @Test
    public static void testCompareByteArrays(){
        int res = SVFastqUtils.compareByteArrays(new byte[] {'A'}, new byte[] {'C'});
        Assert.assertTrue(res<0);
        res = SVFastqUtils.compareByteArrays(new byte[] {'A'}, new byte[] {'A'});
        Assert.assertTrue(res==0);
        res = SVFastqUtils.compareByteArrays(new byte[] {'A', 'C'}, new byte[] {'A', 'C', 'T'});
        Assert.assertTrue(res<0);
        res = SVFastqUtils.compareByteArrays(new byte[] {'A', 'T'}, new byte[] {'A', 'C', 'T'});
        Assert.assertTrue(res>0);
    }

    @Test
    public static void testSortFastqRecords(){
        final List<byte[]> input = Arrays.asList(new byte[] {'A', 'C', 'T'}, new byte[] {'A', 'C'});
        SVFastqUtils.sortFastqRecords(input);
        Assert.assertEquals(input.size(), 2);
        Assert.assertArrayEquals(input.get(0), new byte[] {'A', 'C'});
        Assert.assertArrayEquals(input.get(1), new byte[] {'A', 'C', 'T'});

    }

    @Test
    public static void testReadToFastqRecord(){

        // normal mapped
        byte[] converted = SVFastqUtils.readToFastqRecord(testMappedRead, true);
        List<Integer> breakIndices = new ArrayList<>(4); // guess
        for(int i=0; i<converted.length; ++i){
            if(converted[i]=='\n') breakIndices.add(i);
        }
        Assert.assertTrue(breakIndices.size()==4);
        String customReadName = "@" + testMappedRead.getName() + "/1 mapping=" + testMappedRead.getContig() + ":" + testMappedRead.getAssignedStart() + ";" + testMappedRead.getCigar();
        Assert.assertEquals(new String(Arrays.copyOfRange(converted, 0, breakIndices.get(0))), customReadName);
        Assert.assertArrayEquals(Arrays.copyOfRange(converted, breakIndices.get(0)+1, breakIndices.get(1)), testMappedRead.getBases());
        Assert.assertEquals(new String(Arrays.copyOfRange(converted, breakIndices.get(2)+1, converted.length-1)), ReadUtils.getBaseQualityString(testMappedRead));

        // unmapped
        converted = SVFastqUtils.readToFastqRecord(testUnmappedRead, true);
        breakIndices.clear();
        for(int i=0; i<converted.length; ++i){
            if(converted[i]=='\n') breakIndices.add(i);
        }
        Assert.assertTrue(breakIndices.size()==4);
        customReadName = "@" + testUnmappedRead.getName() + "/1 mapping=unmapped";
        Assert.assertEquals(new String(Arrays.copyOfRange(converted, 0, breakIndices.get(0))), customReadName);
        Assert.assertArrayEquals(Arrays.copyOfRange(converted, breakIndices.get(0)+1, breakIndices.get(1)), testUnmappedRead.getBases());
        Assert.assertEquals(new String(Arrays.copyOfRange(converted, breakIndices.get(2)+1, converted.length-1)), ReadUtils.getBaseQualityString(testUnmappedRead));

        // SA
        converted = SVFastqUtils.readToFastqRecord(testMultiMappedNonprimaryRead, true);
        breakIndices.clear();
        for(int i=0; i<converted.length; ++i){
            if(converted[i]=='\n') breakIndices.add(i);
        }
        Assert.assertTrue(breakIndices.size()==4);
        customReadName = "@" + testMultiMappedNonprimaryRead.getName() + "/1 mapping=" + testMultiMappedNonprimaryRead.getContig() + ":" + testMultiMappedNonprimaryRead.getAssignedStart() + ";" + testMultiMappedNonprimaryRead.getCigar();
        Assert.assertEquals(new String(Arrays.copyOfRange(converted, 0, breakIndices.get(0))), customReadName);
        Assert.assertArrayEquals(Arrays.copyOfRange(converted, breakIndices.get(0)+1, breakIndices.get(1)), testMultiMappedNonprimaryRead.getBases());
        Assert.assertEquals(new String(Arrays.copyOfRange(converted, breakIndices.get(2)+1, converted.length-1)), ReadUtils.getBaseQualityString(testMultiMappedNonprimaryRead));
    }

    @Test
    public static void testConvertFASTQtoGATKRead(){
        final GATKRead read = SVFastqUtils.convertToRead(testMappedReadFastq);
        Assert.assertEquals(read.getName(), testMappedRead.getName() + "/" + (testMappedRead.isFirstOfPair()?"1":"2"));
        Assert.assertArrayEquals(read.getBases(), testMappedRead.getBases());
        Assert.assertEquals(ReadUtils.getBaseQualityString(read), ReadUtils.getBaseQualityString(testMappedRead));
        Assert.assertEquals((int)read.getAttributeAsInteger("RC"), 1);
    }

    @Test
    public static void testWriteFastqFile(){
        final File tempDir =createTempDir("whatever"); tempDir.deleteOnExit();
        final String out = tempDir.toString() + "/test.fastq";
        SVFastqUtils.writeFastqFile(out, PipelineOptionsFactory.create(), Collections.singletonList(SVFastqUtils.readToFastqRecord(testMappedRead, true)));
        final File testFile = new File(tempDir, "test.fastq");
        Assert.assertTrue(testFile.exists());
        Assert.assertTrue(testFile.length()!=0);
    }

    @Test
    public static void testExtractFASTQContents(){
        final List<FastqRecord> list = SVFastqUtils.extractFASTQContents(new String(SVFastqUtils.readToFastqRecord(testMappedRead, true)));
        Assert.assertEquals(list.size(), 1);
        Assert.assertEquals(list.get(0), testMappedReadFastq);
    }

    @Test
    public static void testLoadFASTQFiles(){
        final File tempFile = createTempFile("testFastq", "fastq"); tempFile.deleteOnExit();
        SVFastqUtils.writeFastqFile(tempFile.toString(), PipelineOptionsFactory.create(), Collections.singletonList(SVFastqUtils.readToFastqRecord(testMappedRead, true)));
        final Map<String, String> singleEntryMap = SVFastqUtils.loadFASTQFiles(ctx, tempFile.toString()).collectAsMap();
        Assert.assertEquals(singleEntryMap.size(), 1);
        Assert.assertEquals(singleEntryMap.keySet().iterator().next(), "file:"+tempFile.toString());
        Assert.assertEquals(SVFastqUtils.extractFASTQContents(singleEntryMap.get("file:"+tempFile.toString())).get(0), testMappedReadFastq);
    }
}
