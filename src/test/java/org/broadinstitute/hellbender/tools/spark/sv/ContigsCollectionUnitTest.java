package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.collect.Iterators;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.junit.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.Map;

/**
 * Mostly testing formatting issues.
 */
public final class ContigsCollectionUnitTest extends BaseTest{

    private static JavaSparkContext ctx;

    private static final File testFASTAFile = new File("src/test/resources/org/broadinstitute/hellbender/tools/spark/sv/RunSGAViaProcessBuilderOnSpark/assembly4.pp.ec.filter.pass.merged.rmdup-contigs.fa");

    @BeforeClass
    private static void setupSparkAndTestFile(){
        SparkContextFactory.enableTestSparkContext();
        ctx = SparkContextFactory.getTestSparkContext(Collections.emptyMap());
    }

    @AfterClass
    private static void closeSpark(){
        SparkContextFactory.stopSparkContext(ctx);
    }

    @Test
    public static void testToListOfStrings() throws IOException{

        Assert.assertNull(new ContigsCollection(1, null).toListOfStrings());

        final ContigsCollection collection = new ContigsCollection(1, FileUtils.readLines(testFASTAFile, "UTF-8"));

        Assert.assertEquals(collection.toListOfStrings().size(), Files.lines(testFASTAFile.toPath()).count());
    }

    @Test
    public static void testToPackedFasta() throws IOException{

        final ContigsCollection collection = new ContigsCollection(1, FileUtils.readLines(testFASTAFile, "UTF-8"));

        Assert.assertEquals(StringUtils.countMatches(collection.toPackedFasta(), "|"), Files.lines(testFASTAFile.toPath()).count()-1);
    }

    @Test
    public void testSplitAssemblyLine(){
        Assert.assertFalse(ContigsCollection.splitAssemblyLine("      ").iterator().hasNext());
        final String assemblyLine = ">contig-13 553 0\tATGTATCCCAGAACTTCCAATTAAAAAGAAATAAAAATTTAAAAATTAAAAAGTAAAGCATTTGTCCAGAGTTGTCTGAGGAATTTAGAATACTTACAATGAAATTTTCTCGGCCAGGCGCGGTGGCTCATGCTTGTAATCCCAGCACTTTGGGAGGCCGAGGTGGGCGGATCACAAGGTCAGGAGATCGAGACCAAGATCGAGACCATCCTGGCGAACATGGTGAAACCCCATCTCTACTAAAAATACAAAAAATTAACCAGGCATGGAGGCAGGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCATGAACCCAGGAGGCAGAGCTTGCAGTGAGCCGAGAAGGCGCCACTGCAGTCCAGCCTGGGCGACAGAGCGACACTCCATCTCAAAAAAAAAAAAAAAAGAAAAGAAATTTTCTCTAACACAAATCTATAATATTAAAATCTATAATTGTATCTATGTATTATAAATAGGTTATAAAAATATAAGATCTATATATGCACATATATATAAAAAACATAGAGATACACACATCT";
        final Iterable<Tuple2<String, String>> result = ContigsCollection.splitAssemblyLine(assemblyLine);
        Assert.assertEquals(Iterators.size(result.iterator()),1);
        Assert.assertEquals(result.iterator().next()._1(), ">contig-13 553 0");
        Assert.assertEquals(result.iterator().next()._2(), "ATGTATCCCAGAACTTCCAATTAAAAAGAAATAAAAATTTAAAAATTAAAAAGTAAAGCATTTGTCCAGAGTTGTCTGAGGAATTTAGAATACTTACAATGAAATTTTCTCGGCCAGGCGCGGTGGCTCATGCTTGTAATCCCAGCACTTTGGGAGGCCGAGGTGGGCGGATCACAAGGTCAGGAGATCGAGACCAAGATCGAGACCATCCTGGCGAACATGGTGAAACCCCATCTCTACTAAAAATACAAAAAATTAACCAGGCATGGAGGCAGGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCATGAACCCAGGAGGCAGAGCTTGCAGTGAGCCGAGAAGGCGCCACTGCAGTCCAGCCTGGGCGACAGAGCGACACTCCATCTCAAAAAAAAAAAAAAAAGAAAAGAAATTTTCTCTAACACAAATCTATAATATTAAAATCTATAATTGTATCTATGTATTATAAATAGGTTATAAAAATATAAGATCTATATATGCACATATATATAAAAAACATAGAGATACACACATCT");
    }

    @Test
    public void testLoadContigsCollectionKeyedByAssemblyId() throws IOException{

        final ContigsCollection collection = new ContigsCollection(1, FileUtils.readLines(testFASTAFile, "UTF-8"));
        final File tempFile = createTempFile("whatever", ""); tempFile.deleteOnExit();
        FileUtils.writeStringToFile(tempFile, "1\t"+collection.toPackedFasta());

        final Map<String, ContigsCollection> singleEntryMap = ContigsCollection.loadContigsCollectionKeyedByAssemblyId(ctx, "file://"+tempFile.getAbsolutePath()).collectAsMap();
        Assert.assertEquals(singleEntryMap.keySet().size(), 1);
        Assert.assertEquals(singleEntryMap.values().iterator().next().getContents().size(), 15);
    }

    // below are moved here from AlignAssembledContigsSparkTest.java
    @Test
    public void testFromPackedFasta() throws IOException{
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

        Assert.assertEquals(ContigsCollection.fromPackedFasta(1, packedFasta).getContents().size(), 3);
    }

    @Test
    public void testParseAlignedAssembledContigLine() throws Exception {
        final String line = "100\t>contig-0 2498 0\t1\t7043012\t7044152\t+\t1141M1357S\t60\t1\t1141\t1";
        final AlignmentRegion region1 = ContigsCollection.parseAlignedAssembledContigLine(line);
        Assert.assertEquals(region1.referenceInterval, new SimpleInterval("1", 7043012, 7044152));
        Assert.assertTrue(region1.forwardStrand);
        Assert.assertEquals(region1.forwardStrandCigar.toString(), "1141M1357S");
        Assert.assertEquals(region1.mapqual, 60);
        Assert.assertEquals(region1.startInAssembledContig, 1);
        Assert.assertEquals(region1.endInAssembledContig, 1141);
        Assert.assertEquals(region1.mismatches, 1);

        final String line2 = "100\tcontig-0\t1\t7044151\t7045305\t+\t1343S1155M\t60\t1344\t2498\t3";
        final AlignmentRegion region2 = ContigsCollection.parseAlignedAssembledContigLine(line2);
        Assert.assertEquals(region2.referenceInterval, new SimpleInterval("1", 7044151, 7045305));
        Assert.assertTrue(region2.forwardStrand);
        Assert.assertEquals(region2.forwardStrandCigar.toString(), "1343S1155M");
        Assert.assertEquals(region2.mapqual, 60);
        Assert.assertEquals(region2.startInAssembledContig, 1344);
        Assert.assertEquals(region2.endInAssembledContig, 2498);
        Assert.assertEquals(region2.mismatches, 3);
    }
}
