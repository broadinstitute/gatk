package org.broadinstitute.hellbender.tools.spark.sv;

import com.github.lindenb.jbwa.jni.AlnRgn;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class AlignmentRegionTest {

    @Test
    public void testNormalizeStrandOnCreation() throws Exception {
        final AlignmentRegion ar1p1 = new AlignmentRegion("1", "1", new AlnRgn("1", 96, (byte) '+', "4M4S", 60, 0, 0));
        Assert.assertEquals(ar1p1.referenceInterval, new SimpleInterval("1", 97, 100));
        Assert.assertEquals(ar1p1.startInAssembledContig, 1);
        Assert.assertEquals(ar1p1.endInAssembledContig, 4);

        final AlignmentRegion ar1p2 = new AlignmentRegion("1", "1", new AlnRgn("1", 196, (byte) '-', "4M4H", 60, 0, 0));
        Assert.assertEquals(ar1p2.referenceInterval, new SimpleInterval("1", 197, 200));
        Assert.assertEquals(ar1p2.startInAssembledContig, 5);
        Assert.assertEquals(ar1p2.endInAssembledContig, 8);

        final AlignmentRegion ar2p1 = new AlignmentRegion("1", "1", new AlnRgn("1", 96, (byte) '-', "4M4S", 60, 0, 0));
        Assert.assertEquals(ar2p1.referenceInterval, new SimpleInterval("1", 97, 100));
        Assert.assertEquals(ar2p1.startInAssembledContig, 5);
        Assert.assertEquals(ar2p1.endInAssembledContig, 8);

        final AlignmentRegion ar2p2 = new AlignmentRegion("1", "1", new AlnRgn("1", 196, (byte) '+', "4M4H", 60, 0, 0));
        Assert.assertEquals(ar2p2.referenceInterval, new SimpleInterval("1", 197, 200));
        Assert.assertEquals(ar2p2.startInAssembledContig, 1);
        Assert.assertEquals(ar2p2.endInAssembledContig, 4);

    }

    @Test
    public void testAlignmentRegionOverlap() throws Exception {

        //overlap by 1
        final AlignmentRegion ar1 = new AlignmentRegion("1","1", TextCigarCodec.decode("5M5H"),true,new SimpleInterval("1",1,5),60,1,5,0);
        final AlignmentRegion ar2 = new AlignmentRegion("1","1", TextCigarCodec.decode("5S5M"),true,new SimpleInterval("1",10,16),60,5,10,0);
        Assert.assertEquals(ar1.overlapOnContig(ar2), 1);

        // don't overlap
        final AlignmentRegion ar3 = new AlignmentRegion("1","1", TextCigarCodec.decode("5M5H"),true,new SimpleInterval("1",1,5),60,1,5,0);
        final AlignmentRegion ar4 = new AlignmentRegion("1","1", TextCigarCodec.decode("5S5M"),true,new SimpleInterval("1",11,16),60,6,10,0);
        Assert.assertEquals(ar3.overlapOnContig(ar4), 0);

    }

}