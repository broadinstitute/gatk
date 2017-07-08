package org.broadinstitute.hellbender.tools.spark.sv.sga;

import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.MiniClusterUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Map;

public class ContigCollectionTest extends BaseTest {

    static final String packedFastaWithLength = ">contig-4 521 0|ATTTGTGAGCGCTTTGAGGCCTTTGGAGAAAAAGGAAATGTATTCACAGAAAAACTTGAAAAAAAGCTTCTGGTAAACTGT" +
            "TTTGTAATGTGTACAATCATTTCACAGAGATCAGTGTTTCTTTTCATTGAGCAGCTTGGAAACTCTATTGTTGTAGAATCTGCAAACGGATATTTTTCAGTG" +
            "CTTTGATGCCTGTGTTAAAAAAGGAAATATCCTCACATAAAAAATGACAGAAAATTTATGAGAAACTTCTTTGTGATGTGTGCATTTATGCCACAGAATTGA" +
            "ACCATTTTATGATTGAGCAGTTTGGAAACAGTCTTTTTGTGGAATCTAAAAAGAGATATTTATGAGCGCATTGAGGCCTACAGTAAAAAAGGAAATATCTTC" +
            "ACATAAAAACTAGCAAGGAGCTTTCTAAGAAACTGCTTTGTGATGCGTGAATTCATCTCACAGAGGTAAATGTTTCTTTGCATTGAACAGTGGAAACTCTGT" +
            "TCTTGTAGAATCTGCAAAGTGATATTTGTGAG|>contig-11 164 0|TCACAGAATTGAACGACCTCTTTGTTTGAGCAGTTTGGGAACAGCCTTTTTG" +
            "TAGATTCTGCAAAGGGATATTTGTAAGCCCTTTGAGGACTATGGTGAAAACGTAAATATCTTCACATAACTAGACAGAAGGTTGCTGAAAAGCTGCTTTGTG" +
            "ATGTGTGATT|>contig-0 207 0|GCATTGAACAGTGAAAACTCTGTTCTTGTAGAATCTGCAAAGTGATATTTGTGAGTGTTTTGAGGCCTATGGTGA" +
            "AAAAGGAAATATCTTCAGAAAAACTAGACAGAAGCTTTCTGAGAATATTCTTTGTGATATATGCATTCATCTCACAGAATTGAACGACCTCTTTGTTTGAGC" +
            "AGTTTGGGAACAGCCTTTTTGTAGATTCTG";
    static final String packedFastaWithoutLength = ">contig-104|ATTTGTGAGCGCTTTGAGGCCTTTGGAGAAAAAGGAAATGTATTCACAGAAAAACTTGAAAAAAAGCTTCTGGTAAACTGT" +
            "TTTGTAATGTGTACAATCATTTCACAGAGATCAGTGTTTCTTTTCATTGAGCAGCTTGGAAACTCTATTGTTGTAGAATCTGCAAACGGATATTTTTCAGTG" +
            "CTTTGATGCCTGTGTTAAAAAAGGAAATATCCTCACATAAAAAATGACAGAAAATTTATGAGAAACTTCTTTGTGATGTGTGCATTTATGCCACAGAATTGA" +
            "ACCATTTTATGATTGAGCAGTTTGGAAACAGTCTTTTTGTGGAATCTAAAAAGAGATATTTATGAGCGCATTGAGGCCTACAGTAAAAAAGGAAATATCTTC" +
            "ACATAAAAACTAGCAAGGAGCTTTCTAAGAAACTGCTTTGTGATGCGTGAATTCATCTCACAGAGGTAAATGTTTCTTTGCATTGAACAGTGGAAACTCTGT" +
            "TCTTGTAGAATCTGCAAAGTGATATTTGTGAG|>contig-111|TCACAGAATTGAACGACCTCTTTGTTTGAGCAGTTTGGGAACAGCCTTTTTG" +
            "TAGATTCTGCAAAGGGATATTTGTAAGCCCTTTGAGGACTATGGTGAAAACGTAAATATCTTCACATAACTAGACAGAAGGTTGCTGAAAAGCTGCTTTGTG" +
            "ATGTGTGATT|>contig-100|GCATTGAACAGTGAAAACTCTGTTCTTGTAGAATCTGCAAAGTGATATTTGTGAGTGTTTTGAGGCCTATGGTGA" +
            "AAAAGGAAATATCTTCAGAAAAACTAGACAGAAGCTTTCTGAGAATATTCTTTGTGATATATGCATTCATCTCACAGAATTGAACGACCTCTTTGTTTGAGC" +
            "AGTTTGGGAACAGCCTTTTTGTAGATTCTG";

    public static final String[][] contigsCollectionPackedFastaAndWithAsmId = {
            new String[]{packedFastaWithLength, "1\t" + packedFastaWithLength},
            new String[]{packedFastaWithoutLength, "2\t" + packedFastaWithoutLength}
    };

    @DataProvider(name = "PackedFasta")
    private Object[][] createInputsAndExpectedResults() {
        return contigsCollectionPackedFastaAndWithAsmId;
    }

    @Test(dataProvider = "PackedFasta", groups = "sv")
    public void testFromPackedFasta(final String packedFasta, final String packedFastaAndAsmId){

        ContigsCollection contigsCollection = ContigsCollection.fromPackedFasta(packedFasta);

        Assert.assertEquals(contigsCollection.getContents().size(), 3);
    }

    @Test(dataProvider = "PackedFasta", groups = "sv")
    public void testLoadContigsCollectionKeyedByAssemblyId(final String packedFasta, final String packedFastaAndAsmId) throws Exception {
        MiniClusterUtils.runOnIsolatedMiniCluster(cluster -> {
            //use the HDFS on the mini cluster
            final Path workingDirectory = MiniClusterUtils.getWorkingDir(cluster);
            final Path tempPath = new Path(workingDirectory, "testLocalAssemblyFormatRead");
            final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

            final FileSystem fs = tempPath.getFileSystem(ctx.hadoopConfiguration());
            final FSDataOutputStream fsOutStream = fs.create(tempPath);

            fsOutStream.writeBytes(packedFastaAndAsmId);
            fsOutStream.writeBytes("\n");
            fsOutStream.close();
            fs.deleteOnExit(tempPath);

            Assert.assertTrue(SparkUtils.pathExists(ctx, tempPath));
            final Map<String, ContigsCollection> assemblyIdAndContigs = ContigsCollection.loadContigsCollectionKeyedByAssemblyId(ctx, tempPath.toString()).collectAsMap();
            Assert.assertEquals(assemblyIdAndContigs.size(), 1);
            Assert.assertEquals(assemblyIdAndContigs.keySet().size(), 1);
            Assert.assertEquals(assemblyIdAndContigs.values().size(), 1);
            Assert.assertEquals(assemblyIdAndContigs.values().iterator().next().getContents().size(), 3);
        });
    }
}