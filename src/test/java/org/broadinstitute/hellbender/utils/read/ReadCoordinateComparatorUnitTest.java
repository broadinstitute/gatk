package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public final class ReadCoordinateComparatorUnitTest extends BaseTest{
    @Test
    public void testEqual() throws Exception {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        ReadCoordinateComparator comp= new ReadCoordinateComparator(header);
        final GATKRead r1 = ArtificialReadUtils.createArtificialRead(header, "10M");
        final GATKRead r2 = ArtificialReadUtils.createArtificialRead(header, "10M");
        Assert.assertEquals(comp.compare(r1, r2), 0);
        Assert.assertEquals(comp.compare(r2, r1), 0);
    }

    @Test
    public void testReverse() throws Exception {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        ReadCoordinateComparator comp= new ReadCoordinateComparator(header);
        final GATKRead r1 = ArtificialReadUtils.createArtificialRead(header, "10M");
        r1.setIsReverseStrand(true);
        final GATKRead r2 = ArtificialReadUtils.createArtificialRead(header, "10M");
        Assert.assertEquals(comp.compare(r1, r2), 1);
        Assert.assertEquals(comp.compare(r2, r1), -1);
    }

    @Test
    public void testPaired() throws Exception {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        ReadCoordinateComparator comp= new ReadCoordinateComparator(header);
        final GATKRead r1 = ArtificialReadUtils.createArtificialRead(header, "10M");
        r1.setIsPaired(true);
        r1.setMatePosition(r1.getContig(), r1.getStart() + 10);

        final GATKRead r2 = ArtificialReadUtils.createArtificialRead(header, "10M");
        r2.setIsPaired(true);
        r2.setMatePosition(r2.getContig(), r2.getStart() + 12);
        Assert.assertEquals(comp.compare(r1, r2), -1);
        Assert.assertEquals(comp.compare(r2, r1), 1);
    }

    @Test
    public void testPairedVsUnpaired() throws Exception {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        ReadCoordinateComparator comp= new ReadCoordinateComparator(header);
        final GATKRead r1 = ArtificialReadUtils.createArtificialRead(header, "10M");
        r1.setIsPaired(true);
        r1.setMatePosition(r1.getContig(), r1.getStart() + 10);

        final GATKRead r2 = ArtificialReadUtils.createArtificialRead(header, "10M");
        Assert.assertEquals(comp.compare(r1, r2), 1);
        Assert.assertEquals(comp.compare(r2, r1), -1);
    }

    @Test
    public void testMappedVsUnmapped() throws Exception {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        ReadCoordinateComparator comp= new ReadCoordinateComparator(header);
        final GATKRead r1 = ArtificialReadUtils.createArtificialRead(header, "10M");

        final GATKRead r2 = ArtificialReadUtils.createArtificialRead(header, "10M");
        r2.setIsUnmapped();
        Assert.assertEquals(comp.compare(r1, r2), -1);
        Assert.assertEquals(comp.compare(r2, r1), 1);
    }

    @Test
    public void testBothUnmapped() throws Exception {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        ReadCoordinateComparator comp= new ReadCoordinateComparator(header);
        final GATKRead r1 = ArtificialReadUtils.createArtificialRead(header, "10M");
        r1.setIsUnmapped();

        final GATKRead r2 = ArtificialReadUtils.createArtificialRead(header, "10M");
        r2.setIsUnmapped();
        Assert.assertEquals(comp.compare(r1, r2), 0);
        Assert.assertEquals(comp.compare(r2, r1), 0);
    }

    @Test
    public void testPosition() throws Exception {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        ReadCoordinateComparator comp= new ReadCoordinateComparator(header);
        final GATKRead r1 = ArtificialReadUtils.createArtificialRead(header, "10M");
        r1.setPosition(r1.getContig(), r1.getStart() + 10);
        final GATKRead r2 = ArtificialReadUtils.createArtificialRead(header, "10M");
        Assert.assertEquals(comp.compare(r1, r2), 1);
        Assert.assertEquals(comp.compare(r2, r1), -1);
    }

    @Test
    public void testContigs() throws Exception {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        ReadCoordinateComparator comp= new ReadCoordinateComparator(header);
        final GATKRead r1 = ArtificialReadUtils.createArtificialRead(header, "10M");
        final GATKRead r2 = ArtificialReadUtils.createArtificialRead(header, "10M");
        r2.setPosition(header.getSequence(1).getSequenceName(), r2.getStart());
        Assert.assertEquals(comp.compare(r1, r2), -1);
        Assert.assertEquals(comp.compare(r2, r1), 1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullHeader() throws Exception {
        new ReadCoordinateComparator(null);
    }
}
