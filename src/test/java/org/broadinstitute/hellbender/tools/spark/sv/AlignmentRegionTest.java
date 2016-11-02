package org.broadinstitute.hellbender.tools.spark.sv;

import com.github.lindenb.jbwa.jni.AlnRgn;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

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

    @Test
    public void testGetTotalHardClipping(){
        //"10H10S10D10M10I10S10H"
        final Cigar cigar = new Cigar(Arrays.asList(new CigarElement(10, CigarOperator.H), new CigarElement(10, CigarOperator.S),
                new CigarElement(10, CigarOperator.D), new CigarElement(10, CigarOperator.M),
                new CigarElement(10, CigarOperator.I), new CigarElement(10, CigarOperator.S),
                new CigarElement(10, CigarOperator.H)));
        Assert.assertEquals(AlignmentRegion.getTotalHardClipping(cigar), 20);
    }

    @Test
    public void testStartOfAlignmentInContig(){
        //"10H10S10D10M10I10S10H"
        final Cigar cigar = new Cigar(Arrays.asList(new CigarElement(10, CigarOperator.H), new CigarElement(10, CigarOperator.S),
                new CigarElement(10, CigarOperator.D), new CigarElement(10, CigarOperator.M),
                new CigarElement(10, CigarOperator.I), new CigarElement(10, CigarOperator.S),
                new CigarElement(10, CigarOperator.H)));
        Assert.assertEquals(AlignmentRegion.startOfAlignmentInContig(cigar), 21);
    }

    @Test
    public void testEndOfAlignmentInContig(){
        //"10H10S10D10M10I10S10H"
        final Cigar cigar = new Cigar(Arrays.asList(new CigarElement(10, CigarOperator.H), new CigarElement(10, CigarOperator.S),
                new CigarElement(10, CigarOperator.D), new CigarElement(10, CigarOperator.M),
                new CigarElement(10, CigarOperator.I), new CigarElement(10, CigarOperator.S),
                new CigarElement(10, CigarOperator.H)));
        Assert.assertEquals(AlignmentRegion.endOfAlignmentInContig(70, cigar), 50);
    }

    @Test
    public void testGetNumClippedBases(){
        //"10H10S10D10M10I10S10H"
        final Cigar cigar = new Cigar(Arrays.asList(new CigarElement(10, CigarOperator.H), new CigarElement(10, CigarOperator.S),
                new CigarElement(10, CigarOperator.D), new CigarElement(10, CigarOperator.M),
                new CigarElement(10, CigarOperator.I), new CigarElement(10, CigarOperator.S),
                new CigarElement(10, CigarOperator.H)));
        Assert.assertEquals(AlignmentRegion.getNumClippedBases(true, cigar), 20);
        Assert.assertEquals(AlignmentRegion.getNumClippedBases(false, cigar), 20);
    }
}