package org.broadinstitute.hellbender.tools.walkers.indels;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMValidationError;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLoggerInterface;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GATKReadWriter;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;


public class ConstrainedMateFixingManagerUnitTest extends BaseTest {

    private static SAMFileHeader header;

    @BeforeClass
    public void beforeClass() {
        header = ArtificialReadUtils.createArtificialSamHeader(3, 1, 10000);
    }

    @Test
    public void testCanMoveReads() {
        final int maxInsertSize = 1000;
        final ConstrainedMateFixingManager manager = new ConstrainedMateFixingManager(null, header, maxInsertSize, 1000, 1000);
        final GATKRead readChr0 = ArtificialReadUtils.createArtificialRead(header, "foo", 0, 1, new byte[]{'A'}, new byte[]{'!'}, "1M");
        final GATKRead readChr1 = ArtificialReadUtils.createArtificialRead(header, "foo", 1, 1, 1);

        // before adding/flushig, returns always true
        Assert.assertTrue(manager.canMoveReads(readChr0));
        manager.addRead(readChr0.deepCopy(), false);
        Assert.assertTrue(manager.canMoveReads(readChr1));

        // force to flush
        manager.addRead(readChr1.deepCopy(), false);
        Assert.assertFalse(manager.canMoveReads(readChr0));
        Assert.assertTrue(manager.canMoveReads(readChr1));

        // shift the coordinate to be maxInsertSize apart
        readChr0.setPosition(readChr0.getContig(), readChr0.getStart() + maxInsertSize);
        Assert.assertFalse(manager.canMoveReads(readChr0));
        // shift the coordinate to be maxInsertSize + 1 apart
        readChr0.setPosition(readChr0.getContig(), readChr0.getStart() + 1);
        Assert.assertTrue(manager.canMoveReads(readChr0));
    }

    @Test
    public void testSecondaryAlignmentsDoNotInterfere() {
        final List<GATKRead> properReads = ArtificialReadUtils.createPair(header, "foo", 1, 10, 30, true, false);
        final SAMRecord read1 = properReads.get(0).convertToSAMRecord(header);
        read1.setAlignmentStart(8); // move the read
        read1.setFlags(99);   // first in proper pair, mate negative strand

        final SAMRecord read2Primary = properReads.get(1).convertToSAMRecord(header);
        read2Primary.setFlags(147);   // second in pair, mate unmapped, not primary alignment

        Assert.assertEquals(read1.getInferredInsertSize(), 21);

        final SAMRecord read2NonPrimary = read2Primary.deepCopy();
        read2NonPrimary.setFlags(393);   // second in proper pair, on reverse strand

        final ConstrainedMateFixingManager manager = new ConstrainedMateFixingManager(null, header, 1000, 1000, 1000);
        manager.addRead(new SAMRecordToGATKReadAdapter(read1), true, false);
        manager.addRead(new SAMRecordToGATKReadAdapter(read2NonPrimary), false, false);
        manager.addRead(new SAMRecordToGATKReadAdapter(read2Primary), false, false);

        Assert.assertEquals(manager.getNReadsInQueue(), 3);

        for ( final GATKRead gatkRead : manager.getReadsInQueueForTesting() ) {
            final SAMRecord read = gatkRead.convertToSAMRecord(header);
            if ( read.getFirstOfPairFlag() ) {
                Assert.assertEquals(read.getFlags(), 99);
                Assert.assertEquals(read.getInferredInsertSize(), 23);
            } else if ( read.getNotPrimaryAlignmentFlag() ) {
                Assert.assertEquals(read.getFlags(), 393);
                Assert.assertEquals(read.getInferredInsertSize(), -21);
            } else {
                Assert.assertEquals(read.getFlags(), 147);
                Assert.assertEquals(read.getInferredInsertSize(), -23);
            }
        }
    }

    @Test
    public void testSecondaryAlignmentsDoNotCauseAccidentalRemovalOfMate() {
        final List<GATKRead> properReads = ArtificialReadUtils.createPair(header, "foo", 1, 530, 1594, true, false);
        final SAMRecord read1 = properReads.get(0).convertToSAMRecord(header);
        read1.setFlags(99);   // first in proper pair, mate negative strand

        final SAMRecord read2Primary = properReads.get(1).convertToSAMRecord(header);
        read2Primary.setFlags(147);   // second in pair, mate unmapped, not primary alignment
        read2Primary.setAlignmentStart(1596); // move the read

        final SAMRecord read2NonPrimary = read2Primary.deepCopy();
        read2NonPrimary.setReadName("foo");
        read2NonPrimary.setFlags(393);   // second in proper pair, on reverse strand
        read2NonPrimary.setAlignmentStart(451);
        read2NonPrimary.setMateAlignmentStart(451);

        final ConstrainedMateFixingManager manager = new ConstrainedMateFixingManager(null, header, 10000, 200, 10000);
        manager.addRead(new SAMRecordToGATKReadAdapter(read2NonPrimary), false, false);
        manager.addRead(new SAMRecordToGATKReadAdapter(read1), false, false);

        for ( int i = 0; i < ConstrainedMateFixingManager.EMIT_FREQUENCY; i++ )
            manager.addRead(ArtificialReadUtils.createArtificialRead(header, "foo" + i, 0, 1500, 10), false, false);

        Assert.assertTrue(manager.forMateMatching.containsKey("foo"));
    }

    @Test
    public void testSupplementaryAlignmentsDoNotCauseBadMateFixing() {
        final List<GATKRead> properReads = ArtificialReadUtils.createPair(header, "foo", 1, 1000, 2000, true, false);
        final SAMRecord read1 = properReads.get(0).convertToSAMRecord(header);
        read1.setFlags(99);   // first in pair, negative strand

        final SAMRecord read2 = properReads.get(1).convertToSAMRecord(header);
        read2.setFlags(161);   // second in pair, mate negative strand

        final SAMRecord read2Supp = read2.deepCopy();
        read2Supp.setReadName("foo");
        read2Supp.setFlags(2209);   // second in pair, mate negative strand, supplementary
        read2Supp.setAlignmentStart(100);
        read2Supp.setMateAlignmentStart(1000);

        final DummyWriter writer = new DummyWriter();
        final ConstrainedMateFixingManager manager = new ConstrainedMateFixingManager(writer, header, 10000, 200, 10000);
        manager.addRead(new SAMRecordToGATKReadAdapter(read2Supp), false, false);
        manager.addRead(new SAMRecordToGATKReadAdapter(read1), false, false);
        manager.addRead(new SAMRecordToGATKReadAdapter(read2), false, false);
        manager.close(); // "write" the reads to our dummy writer

        // check to make sure that none of the mate locations were changed, which is the problem brought to us by a user
        for ( final GATKRead gatkRead : writer.reads ) {
            final SAMRecord read = gatkRead.convertToSAMRecord(header);
            final int start = read.getAlignmentStart();
            switch (start) {
                case 100:
                    Assert.assertEquals(read.getMateAlignmentStart(), 1000);
                    break;
                case 1000:
                    Assert.assertEquals(read.getMateAlignmentStart(), 2000);
                    break;
                case 2000:
                    Assert.assertEquals(read.getMateAlignmentStart(), 1000);
                    break;
                default:
                    Assert.assertTrue(false, "We saw a read located at the wrong position");
            }
        }
    }

    private class DummyWriter implements GATKReadWriter {

        public List<GATKRead> reads;

        public DummyWriter() { reads = new ArrayList<>(10); }

        public void addRead(final GATKRead alignment) { reads.add(alignment);}

        public void close() {}
    }
}