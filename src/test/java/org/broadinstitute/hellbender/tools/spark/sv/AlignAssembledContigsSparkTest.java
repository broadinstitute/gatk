package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.commons.collections4.IterableUtils;
import org.broadinstitute.hellbender.tools.spark.sv.ContigAligner.AlignmentRegion;
import org.broadinstitute.hellbender.tools.spark.sv.ContigAligner.AssembledBreakpoint;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.util.List;

public class AlignAssembledContigsSparkTest extends BaseTest {

    @Test
    public void testPackedFasta() throws Exception {
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

        RunSGAViaProcessBuilderOnSpark.ContigsCollection contigsCollection = RunSGAViaProcessBuilderOnSpark.ContigsCollection.fromPackedFasta(packedFasta);

        Assert.assertEquals(contigsCollection.getContents().size(), 3);
    }

    @Test
    public void testParseAlignedAssembledContigLine() throws Exception {
        final String line = "100\t>contig-0 2498 0\t1\t7043012\t7044153\t+\t1141M1357S\t60\t1\t1141\t1";
        final AlignmentRegion region1 = AlignAssembledContigsSpark.parseAlignedAssembledContigLine(line);
        Assert.assertEquals(region1.referenceInterval, new SimpleInterval("1", 7043012, 7044153));
        Assert.assertTrue(region1.forwardStrand);
        Assert.assertEquals(region1.cigar.toString(), "1141M1357S");
        Assert.assertEquals(region1.mqual, 60);
        Assert.assertEquals(region1.startInAssembledContig, 1);
        Assert.assertEquals(region1.endInAssembledContig, 1141);
        Assert.assertEquals(region1.mismatches, 1);

        final String line2 = "100\tcontig-0\t1\t7044151\t7045306\t+\t1343S1155M\t60\t1344\t2498\t3";
        final AlignmentRegion region2 = AlignAssembledContigsSpark.parseAlignedAssembledContigLine(line2);
        Assert.assertEquals(region2.referenceInterval, new SimpleInterval("1", 7044151, 7045306));
        Assert.assertTrue(region2.forwardStrand);
        Assert.assertEquals(region2.cigar.toString(), "1343S1155M");
        Assert.assertEquals(region2.mqual, 60);
        Assert.assertEquals(region2.startInAssembledContig, 1344);
        Assert.assertEquals(region2.endInAssembledContig, 2498);
        Assert.assertEquals(region2.mismatches, 3);
    }
}