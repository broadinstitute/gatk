package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public final class ReadCoordinateComparatorUnitTest extends GATKBaseTest {

    /**
     * Tests that the ordering produced by ReadCoordinateComparator matches coordinate file ordering
     * as produced by htsjdk's {@link SAMRecordCoordinateComparator} for a representative selection of reads. Ignores
     * differences in tie-breaking done for reads with the same position -- just asserts that the reads are
     * coordinate-sorted according to htsjdk, including unmapped reads with and without an assigned position.
     */
    @Test
    public void testComparatorOrderingMatchesHtsjdkFileOrdering() throws IOException {
        // This file has unmapped reads that are set to the position of their mates, as well as
        // unmapped reads with no position -- the ordering check in the test below will fail if
        // our ordering of these reads relative to the mapped reads is not consistent with the
        // definition of coordinate sorting as defined in htsjdk.samtools.SAMRecordCoordinateComparator
        final String inputBam = publicTestDir + "org/broadinstitute/hellbender/utils/read/comparator_test_with_unmapped.bam";
        final List<GATKRead> reads = new ArrayList<>();
        SAMFileHeader header = null;

        try ( final ReadsDataSource readsSource = new ReadsDataSource(IOUtils.getPath(inputBam)) ) {
            header = readsSource.getHeader();

            for ( GATKRead read : readsSource ) {
                reads.add(read);
            }
        }

        // Randomize ordering and then re-sort with the ReadCoordinateComparator
        Collections.shuffle(reads);
        Collections.sort(reads, new ReadCoordinateComparator(header));

        // Check ordering using SAMRecordCoordinateComparator.fileOrderCompare() instead of the full comparator,
        // to ignore meaningless differences in tie-breaking rules for reads with the same position.
        final SAMRecordCoordinateComparator samComparator = new SAMRecordCoordinateComparator();
        GATKRead previousRead = null;
        for ( final GATKRead currentRead : reads ) {
            if ( previousRead != null ) {
                Assert.assertTrue(samComparator.fileOrderCompare(previousRead.convertToSAMRecord(header), currentRead.convertToSAMRecord(header)) <= 0,
                                                                 "Reads are out of order: " + previousRead + " and " + currentRead);
            }
            previousRead = currentRead;
        }
    }

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
        ReadCoordinateComparator comp = new ReadCoordinateComparator(header);

        final GATKRead r1 = ArtificialReadUtils.createArtificialRead(header, "10M");
        final GATKRead r2 = ArtificialReadUtils.createArtificialUnmappedRead(header, new byte[]{'A'}, new byte[]{30});
        Assert.assertTrue(! r1.isUnmapped() && r2.isUnmapped());

        Assert.assertEquals(comp.compare(r1, r2), -1);
        Assert.assertEquals(comp.compare(r2, r1), 1);
    }

    @Test
    public void testBothUnmapped() throws Exception {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        ReadCoordinateComparator comp= new ReadCoordinateComparator(header);

        final GATKRead r1 = ArtificialReadUtils.createArtificialUnmappedRead(header, new byte[]{'A'}, new byte[]{30});
        final GATKRead r2 = ArtificialReadUtils.createArtificialUnmappedRead(header, new byte[]{'A'}, new byte[]{30});
        Assert.assertTrue(r1.isUnmapped() && r2.isUnmapped());

        Assert.assertEquals(comp.compare(r1, r2), 0);
        Assert.assertEquals(comp.compare(r2, r1), 0);
    }

    @Test
    public void testUnmappedWithAssignedPositionVsMapped() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        ReadCoordinateComparator comp= new ReadCoordinateComparator(header);

        final GATKRead r1 = ArtificialReadUtils.createArtificialRead(header, "10M");
        final GATKRead r2 = ArtificialReadUtils.createArtificialUnmappedReadWithAssignedPosition(header,r1.getContig(), r1.getStart() - 1, new byte[]{'A'}, new byte[]{30});
        final GATKRead r3 = ArtificialReadUtils.createArtificialUnmappedReadWithAssignedPosition(header,r1.getContig(), r1.getStart() + 1, new byte[]{'A'}, new byte[]{30});
        Assert.assertTrue((! r1.isUnmapped()) && r2.isUnmapped() && r3.isUnmapped());

        Assert.assertEquals(comp.compare(r1, r2), 1);
        Assert.assertEquals(comp.compare(r2, r1), -1);
        Assert.assertEquals(comp.compare(r1, r3), -1);
        Assert.assertEquals(comp.compare(r3, r1), 1);
        Assert.assertEquals(comp.compare(r2, r3), -1);
        Assert.assertEquals(comp.compare(r3, r2), 1);
    }

    @Test
    public void testUnmappedWithAssignedPositionVsUnmapped() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        ReadCoordinateComparator comp = new ReadCoordinateComparator(header);

        final GATKRead r1 = ArtificialReadUtils.createArtificialUnmappedReadWithAssignedPosition(header, "1", 1, new byte[]{'A'}, new byte[]{30});
        final GATKRead r2 = ArtificialReadUtils.createArtificialUnmappedRead(header, new byte[]{'A'}, new byte[]{30});
        Assert.assertTrue(r1.isUnmapped() && r2.isUnmapped());

        Assert.assertEquals(comp.compare(r1, r2), -1);
        Assert.assertEquals(comp.compare(r2, r1), 1);
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
