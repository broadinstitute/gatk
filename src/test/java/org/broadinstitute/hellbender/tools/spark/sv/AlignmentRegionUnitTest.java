package org.broadinstitute.hellbender.tools.spark.sv;

import com.github.lindenb.jbwa.jni.AlnRgn;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

public class AlignmentRegionUnitTest {

    @Test
    public void testNormalizeStrandOnCreation() throws Exception {
        final AlignmentRegion ar1p1 = new AlignmentRegion("1", "1", new AlnRgn("1", 96, (byte) '+', "4M4S", 60, 0, 0));
        Assert.assertEquals(ar1p1.referenceInterval, new SimpleInterval("1", 97, 100));
        Assert.assertEquals(ar1p1.startInAssembledContig, 1);
        Assert.assertEquals(ar1p1.endInAssembledContig, 4);
        Assert.assertTrue(ar1p1.forwardStrand);
        Assert.assertEquals(ar1p1.forwardStrandCigar, TextCigarCodec.decode("4M4S"));
        Assert.assertEquals(ar1p1.assembledContigLength, 8);
        Assert.assertEquals(ar1p1.mismatches, 0);

        final AlignmentRegion ar1p2 = new AlignmentRegion("1", "1", new AlnRgn("1", 196, (byte) '-', "4M4H", 60, 1, 0));
        Assert.assertEquals(ar1p2.referenceInterval, new SimpleInterval("1", 197, 200));
        Assert.assertEquals(ar1p2.startInAssembledContig, 5);
        Assert.assertEquals(ar1p2.endInAssembledContig, 8);
        Assert.assertFalse(ar1p2.forwardStrand);
        Assert.assertEquals(ar1p2.forwardStrandCigar, TextCigarCodec.decode("4H4M"));
        Assert.assertEquals(ar1p2.assembledContigLength, 8);
        Assert.assertEquals(ar1p2.mismatches, 1);

        final AlignmentRegion ar2p1 = new AlignmentRegion("1", "1", new AlnRgn("1", 96, (byte) '-', "4M4S", 60, 0, 0));
        Assert.assertEquals(ar2p1.referenceInterval, new SimpleInterval("1", 97, 100));
        Assert.assertEquals(ar2p1.startInAssembledContig, 5);
        Assert.assertEquals(ar2p1.endInAssembledContig, 8);
        Assert.assertFalse(ar2p1.forwardStrand);
        Assert.assertEquals(ar2p1.forwardStrandCigar, TextCigarCodec.decode("4S4M"));
        Assert.assertEquals(ar2p1.assembledContigLength, 8);
        Assert.assertEquals(ar2p1.mismatches, 0);

        final AlignmentRegion ar2p2 = new AlignmentRegion("1", "1", new AlnRgn("1", 196, (byte) '+', "4M4H", 60, 0, 0));
        Assert.assertEquals(ar2p2.referenceInterval, new SimpleInterval("1", 197, 200));
        Assert.assertEquals(ar2p2.startInAssembledContig, 1);
        Assert.assertEquals(ar2p2.endInAssembledContig, 4);
        Assert.assertTrue(ar2p2.forwardStrand);
        Assert.assertEquals(ar2p2.forwardStrandCigar, TextCigarCodec.decode("4M4H"));
        Assert.assertEquals(ar2p1.assembledContigLength, 8);
        Assert.assertEquals(ar2p2.mismatches, 0);
    }

    @Test
    public void testCtorConcordance() {
        final AlignmentRegion ar1p1 = new AlignmentRegion("1", "1", new AlnRgn("1", 96, (byte) '+', "4M4S", 60, 0, 0));
        final AlignmentRegion ar1p1v2 = new AlignmentRegion("1", "1", TextCigarCodec.decode("4M4S"), true, new SimpleInterval("1", 97, 100), 60, 1, 4, 0);
        Assert.assertEquals(ar1p1, ar1p1v2);

        final AlignmentRegion ar1p2 = new AlignmentRegion("1", "1", new AlnRgn("1", 196, (byte) '-', "4M4H", 60, 0, 0));
        final AlignmentRegion ar1p2v2 = new AlignmentRegion("1", "1", TextCigarCodec.decode("4H4M"), false, new SimpleInterval("1", 197, 200), 60, 5, 8, 0);
        Assert.assertEquals(ar1p2, ar1p2v2);
    }

    @Test
    public void testParseAlignedAssembledContigLine() throws Exception {
        final String line = "100\t>contig-0 2498 0\t1\t7043012\t7044153\t+\t1141M1357S\t60\t1\t1141\t1";
        final AlignmentRegion region1 = AlignmentRegion.parseAlignedAssembledContigLine(line);
        Assert.assertEquals(region1.referenceInterval, new SimpleInterval("1", 7043012, 7044153));
        Assert.assertTrue(region1.forwardStrand);
        Assert.assertEquals(region1.forwardStrandCigar.toString(), "1141M1357S");
        Assert.assertEquals(region1.mapQual, 60);
        Assert.assertEquals(region1.startInAssembledContig, 1);
        Assert.assertEquals(region1.endInAssembledContig, 1141);
        Assert.assertEquals(region1.mismatches, 1);

        final String line2 = "100\tcontig-0\t1\t7044151\t7045306\t+\t1343S1155M\t60\t1344\t2498\t3";
        final AlignmentRegion region2 = AlignmentRegion.parseAlignedAssembledContigLine(line2);
        Assert.assertEquals(region2.referenceInterval, new SimpleInterval("1", 7044151, 7045306));
        Assert.assertTrue(region2.forwardStrand);
        Assert.assertEquals(region2.forwardStrandCigar.toString(), "1343S1155M");
        Assert.assertEquals(region2.mapQual, 60);
        Assert.assertEquals(region2.startInAssembledContig, 1344);
        Assert.assertEquals(region2.endInAssembledContig, 2498);
        Assert.assertEquals(region2.mismatches, 3);
    }

    @Test
    public void testPositionsOnContig(){
        final AlignmentRegion ar1p1 = new AlignmentRegion("1", "1", new AlnRgn("1", 96, (byte) '+', "4M4S", 60, 0, 0));
        Assert.assertEquals(ar1p1.startOfAlignmentInContig(), 1);
        Assert.assertEquals(ar1p1.endOfAlignmentInContig(), 4);

        final AlignmentRegion ar1p2 = new AlignmentRegion("1", "1", new AlnRgn("1", 196, (byte) '-', "4M4H", 60, 0, 0));
        Assert.assertEquals(ar1p2.startOfAlignmentInContig(), 5);
        Assert.assertEquals(ar1p2.endOfAlignmentInContig(), 8);

        final AlignmentRegion ar2p1 = new AlignmentRegion("1", "1", new AlnRgn("1", 96, (byte) '-', "4M4S", 60, 0, 0));
        Assert.assertEquals(ar2p1.startOfAlignmentInContig(), 5);
        Assert.assertEquals(ar2p1.endOfAlignmentInContig(), 8);

        final AlignmentRegion ar2p2 = new AlignmentRegion("1", "1", new AlnRgn("1", 196, (byte) '+', "4M4H", 60, 0, 0));
        Assert.assertEquals(ar2p2.startOfAlignmentInContig(), 1);
        Assert.assertEquals(ar2p2.endOfAlignmentInContig(), 4);
    }
}