package org.broadinstitute.hellbender.tools.walkers.rnaseq;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GATKReadWriter;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class OverhangFixingManagerUnitTest extends GATKBaseTest {

    private SAMFileHeader getHG19Header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(hg19GenomeLocParser.getSequenceDictionary());
        return header;
    }


    @Test
    public void testCleanSplices() {

        final OverhangFixingManager manager = new OverhangFixingManager(getHG19Header(), null, hg19GenomeLocParser, hg19ReferenceReader, 10000, 1, 40, false, true);
        manager.activateWriting();
        final int offset = 10;
        for ( int i = 0; i < OverhangFixingManager.MAX_SPLICES_TO_KEEP + 1; i++ )
            manager.addSplicePosition("1", offset + i, offset + 1 + i);

        final List<OverhangFixingManager.Splice> splices = manager.getSplicesForTesting();

        Assert.assertEquals(splices.size(), (OverhangFixingManager.MAX_SPLICES_TO_KEEP / 2) + 1);

        final int minStartPos = (OverhangFixingManager.MAX_SPLICES_TO_KEEP / 2) + offset;

        for ( final OverhangFixingManager.Splice splice : splices )
            Assert.assertTrue(splice.loc.getStart() >= minStartPos);
    }

    @DataProvider(name = "OverhangTest")
    public Object[][] makeOverhangData() {
        final List<Object[]> tests = new ArrayList<>();
        for ( int leftRead : Arrays.asList(10, 20, 30, 40) ) {
            for ( int rightRead : Arrays.asList(20, 30, 40, 50) ) {
                if ( leftRead >= rightRead )
                    continue;
                for ( int leftSplice : Arrays.asList(10, 20, 30) ) {
                    for ( int rightSplice : Arrays.asList(20, 30, 40) ) {
                        if ( leftSplice >= rightSplice )
                            continue;

                        final GenomeLoc readLoc = hg19GenomeLocParser.createGenomeLoc("1", leftRead, rightRead);
                        final GenomeLoc spliceLoc = hg19GenomeLocParser.createGenomeLoc("1", leftSplice, rightSplice);
                        tests.add(new Object[]{readLoc, spliceLoc});
                    }
                }
            }
        }
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "OverhangTest")
    public void testLeftOverhangs(final GenomeLoc readLoc, final GenomeLoc spliceLoc) {
        final boolean isValidOverhang = readLoc.getStart() <= spliceLoc.getStop() &&
                readLoc.getStop() > spliceLoc.getStop() &&
                readLoc.getStart() > spliceLoc.getStart();
        Assert.assertEquals(OverhangFixingManager.isLeftOverhang(readLoc, spliceLoc), isValidOverhang, readLoc + " vs. " + spliceLoc);
    }

    @Test(dataProvider = "OverhangTest")
    public void testRightOverhangs(final GenomeLoc readLoc, final GenomeLoc spliceLoc) {
        final boolean isValidOverhang = readLoc.getStop() >= spliceLoc.getStart() &&
                readLoc.getStop() < spliceLoc.getStop() &&
                readLoc.getStart() < spliceLoc.getStart();
        Assert.assertEquals(OverhangFixingManager.isRightOverhang(readLoc, spliceLoc), isValidOverhang, readLoc + " vs. " + spliceLoc);
    }

    @DataProvider(name = "MismatchEdgeConditionTest")
    public Object[][] makeMismatchEdgeConditionData() {
        final List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{null, 1, null, 1, 0});
        tests.add(new Object[]{null, 1, null, 1, 100});
        tests.add(new Object[]{new byte[4], 1, null, 1, 3});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MismatchEdgeConditionTest")
    public void testMismatchEdgeCondition(final byte[] read, final int readStart, final byte[] ref, final int refStart, final int overhang) {
        final OverhangFixingManager manager = new OverhangFixingManager(getHG19Header(), null, hg19GenomeLocParser, hg19ReferenceReader, 10000, 1, 40, false, true);
        manager.activateWriting();
        Assert.assertFalse(manager.overhangingBasesMismatch(read, readStart, ((read==null)?0:read.length), ref, refStart, overhang));
    }

    @DataProvider(name = "MismatchTest")
    public Object[][] makeMismatchData() {
        final List<Object[]> tests = new ArrayList<>();

        final byte[] AAAA = new byte[]{(byte)'A', (byte)'A', (byte)'A', (byte)'A'};
        final byte[] AAAC = new byte[]{(byte)'A', (byte)'A', (byte)'A', (byte)'C'};
        final byte[] AAAAAA = new byte[]{(byte)'A', (byte)'A', (byte)'A', (byte)'A', (byte)'A', (byte)'A'};
        final byte[] AAAACA = new byte[]{(byte)'A', (byte)'A', (byte)'A', (byte)'A', (byte)'C', (byte)'A'};
        final byte[] AAAACC = new byte[]{(byte)'A', (byte)'A', (byte)'A', (byte)'A', (byte)'C', (byte)'C'};
        final byte[] AACCAA = new byte[]{(byte)'A', (byte)'A', (byte)'C', (byte)'C', (byte)'A', (byte)'A'};

        tests.add(new Object[]{AAAA, 2, AAAA, 2, 2, false});
        tests.add(new Object[]{AAAA, 2, AAAC, 2, 2, true});
        tests.add(new Object[]{AAAAAA, 3, AAAACA, 3, 3, false});
        tests.add(new Object[]{AAAAAA, 3, AAAACC, 3, 3, true});
        tests.add(new Object[]{AAAAAA, 4, AAAACC, 4, 2, true});
        tests.add(new Object[]{AAAAAA, 2, AAAACC, 2, 3, false});
        tests.add(new Object[]{AAAAAA, 0, AACCAA, 4, 1, false});
        tests.add(new Object[]{AAAAAA, 0, AACCAA, 4, 2, false});
        tests.add(new Object[]{AAAAAA, 3, AACCAA, 0, 2, false});
        tests.add(new Object[]{AAAAAA, 3, AACCAA, 0, 3, false});
        tests.add(new Object[]{AAAAAA, 2, AACCAA, 0, 4, false}); //Testing the length filter

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MismatchTest")
    public void testMismatch(final byte[] read, final int readStart, final byte[] ref, final int refStart, final int overhang, final boolean expected) {
        final OverhangFixingManager manager = new OverhangFixingManager(getHG19Header(), null, hg19GenomeLocParser, hg19ReferenceReader, 10000, 1, 40, false, true);
        manager.activateWriting();
        Assert.assertEquals(manager.overhangingBasesMismatch(read, readStart, ((read==null)?0:read.length), ref, refStart, overhang), expected, new String(read) + " vs. " + new String(ref) + " @" + overhang);
    }

    @DataProvider(name = "unmappedRead")
    public Object[][] makeUnmappedRead() {
        final List<Object[]> tests = new ArrayList<>();
        final GATKRead read = ArtificialReadUtils.createRandomRead(100);
        read.setName("foo");
        read.setCigar("*");
        read.setIsUnmapped();

        tests.add(new Object[]{read});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "unmappedRead")
     public void testUnmappedReadsDoNotFail(final GATKRead read) {
        // try to add unmapped read to the manager
        final OverhangFixingManager manager = new OverhangFixingManager(getHG19Header(), null, null, null, 100, 1, 30, false, true);
        manager.activateWriting();
        manager.addReadGroup(Collections.singletonList(read)); // we just want to make sure that the following call does not fail
    }


   /* This tests that if we peek in the waiting reads (which occurs if there are more than one read added) we don't get
    an NPE from getContig() of an unmapped read.*/
    @Test(dataProvider = "unmappedRead")
    public void testTwoUnmappedReadsDoNotFail(final GATKRead read) {
        // try to add unmapped read to the manager twice
        final OverhangFixingManager manager = new OverhangFixingManager(getHG19Header(), null, null, null, 100, 1, 30, false, true);
        manager.activateWriting();
        manager.addReadGroup(Collections.singletonList(read)); // we need to check when there is an unmapped read waiting to be added (so we need to add two)
        manager.addReadGroup(Collections.singletonList(read)); // we just want to make sure that the following call does not fail
    }

    @Test
    public void testMappingReadMateRepair() {
        final OverhangFixingManager manager = new OverhangFixingManagerAlwaysSplit10000Reads(getHG19Header(), null, hg19GenomeLocParser, hg19ReferenceReader, 10000, 1, 40, false, false);
        GATKRead read1a = ArtificialReadUtils.createArtificialRead(hg19Header, "read1", 1, 10000, new byte[]{(byte)'A', (byte)'A', (byte)'A', (byte)'A', (byte)'A', (byte)'A'}, new byte[]{(byte)'A', (byte)'A', (byte)'A', (byte)'A', (byte)'A', (byte)'A'}, "6M");
        read1a.setMatePosition("1", 10020);
        read1a.setAttribute("MC","6M");
        read1a.setIsFirstOfPair();
        GATKRead read1b = ArtificialReadUtils.createArtificialRead(hg19Header, "read1", 1, 10010, new byte[]{(byte)'C', (byte)'C', (byte)'A', (byte)'A', (byte)'A', (byte)'A'}, new byte[]{(byte)'C', (byte)'C', (byte)'A', (byte)'A', (byte)'A', (byte)'A'}, "6M");
        read1b.setMatePosition("1",10020);
        read1b.setAttribute("MC","6M");
        read1b.setIsFirstOfPair();
        List<GATKRead> readGroup = new LinkedList<GATKRead>();
        readGroup.add(read1a);
        readGroup.add(read1b);

        GATKRead read2 = ArtificialReadUtils.createArtificialRead(hg19Header, "read1", 1, 10020, new byte[]{(byte)'C', (byte)'C', (byte)'A', (byte)'A', (byte)'A', (byte)'A'}, new byte[]{(byte)'C', (byte)'C', (byte)'A', (byte)'A', (byte)'A', (byte)'A'}, "6M");
        read2.setMatePosition("1",10000);
        read2.setAttribute("MC","6M");
        read2.setIsSecondOfPair();

        manager.addReadGroup(readGroup);
        manager.addSplicePosition("1", 1, 2); // this will activate the splitting code
        manager.activateWriting(); // this will flush and populate the repair tree

        Assert.assertFalse(manager.setPredictedMateInformation(read1a));
        Assert.assertEquals(read1a.getAttributeAsString("MC"), "6M");
        Assert.assertEquals(read1a.getMateStart(), 10020);
        Assert.assertFalse(manager.setPredictedMateInformation(read1b));
        Assert.assertEquals(read1b.getAttributeAsString("MC"), "6M");
        Assert.assertEquals(read1b.getMateStart(), 10020);
        Assert.assertTrue(manager.setPredictedMateInformation(read2));
        Assert.assertEquals(read2.getAttributeAsString("MC"), "3S3M");
        Assert.assertEquals(read2.getMateStart(), 10003);
    }

    @Test
    // Test asserting that the failure found in https://github.com/broadinstitute/gatk/issues/6776 has been fixed
    public void testIncidentalMatchFromReportedFailure() {
        GATKRead read1 = ArtificialReadUtils.createArtificialRead(hg19Header, "RR.6123", 1, 3500, 100);
        read1.setIsFirstOfPair();
        GATKRead read2 = ArtificialReadUtils.createArtificialUnmappedRead(hg19Header, new byte[]{(byte)'A'}, new byte[]{(byte)'A'});
        read2.setName("RR.6123135");
        read2.setIsSecondOfPair();

        // Asserting that in the case of reads that end with widely ranged numbers it is no longer possible to incidentally overlap in key generation
        Assert.assertNotEquals( OverhangFixingManager.makeKey(read1.getName(), read1.isFirstOfPair(), read1.getStart()),
                                OverhangFixingManager.makeKey(read2.getName(), read2.isFirstOfPair(), read2.getStart()));
    }

    @Test
    public void testMappingReadMateRepairNoMCTag() {
        final OverhangFixingManager manager = new OverhangFixingManagerAlwaysSplit10000Reads(getHG19Header(), null, hg19GenomeLocParser, hg19ReferenceReader, 10000, 1, 40, false, false);
        GATKRead read1a = ArtificialReadUtils.createArtificialRead(hg19Header, "read1", 1, 10000, "AAAAAA".getBytes(), "AAAAAA".getBytes(), "6M");
        read1a.setMatePosition("1", 10020);
        read1a.setIsFirstOfPair();
        GATKRead read1b = ArtificialReadUtils.createArtificialRead(hg19Header, "read1", 1, 10010, "CCAAAA".getBytes(), "CCAAAA".getBytes(), "6M");
        read1b.setMatePosition("1",10020);
        read1b.setIsFirstOfPair();
        List<GATKRead> readGroup1 = new LinkedList<GATKRead>();
        readGroup1.add(read1a);
        readGroup1.add(read1b);
        GATKRead read2 = ArtificialReadUtils.createArtificialRead(hg19Header, "read1", 1, 10020,  "CCAAAA".getBytes(),  "CCAAAA".getBytes(), "6M");
        read2.setMatePosition("1",10000);
        read2.setIsSecondOfPair();

        // Second Alignment set as secondary to see that it wont be clipped
        GATKRead read3a = ArtificialReadUtils.createArtificialRead(hg19Header, "read2", 1, 10000,"AAAAAA".getBytes(), "AAAAAA".getBytes(), "6M");
        read3a.setMatePosition("1", 10020);
        read3a.setIsSecondaryAlignment(true);
        read3a.setIsSecondOfPair();
        List<GATKRead> readGroup2 = Collections.singletonList(read3a);
        GATKRead read4 = ArtificialReadUtils.createArtificialRead(hg19Header, "read2", 1, 10020, "CCAAAA".getBytes(), "CCAAAA".getBytes(), "6M");
        read4.setMatePosition("1",10000);
        read4.setIsFirstOfPair();

        manager.addReadGroup(readGroup1);
        manager.addReadGroup(readGroup2);
        manager.addSplicePosition("1", 1, 2); // this will activate the splitting code
        manager.activateWriting(); // this will flush and populate the repair tree

        Assert.assertFalse(manager.setPredictedMateInformation(read1a));
        Assert.assertEquals(read1a.getMateStart(), 10020);
        Assert.assertFalse(manager.setPredictedMateInformation(read1b));
        Assert.assertEquals(read1b.getMateStart(), 10020);
        Assert.assertTrue(manager.setPredictedMateInformation(read2));
        Assert.assertEquals(read2.getMateStart(), 10003);

        // Despite having been split, read 3a should not have been marked as changing
        Assert.assertFalse(manager.setPredictedMateInformation(read3a));
        Assert.assertFalse(manager.setPredictedMateInformation(read4));
    }

    private static class OverhangFixingManagerAlwaysSplit10000Reads extends OverhangFixingManager {
        public OverhangFixingManagerAlwaysSplit10000Reads(SAMFileHeader header, GATKReadWriter writer, GenomeLocParser genomeLocParser, ReferenceSequenceFile referenceReader, int maxRecordsInMemory, int maxMismatchesInOverhangs, int maxBasesInOverhangs, boolean doNotFixOverhangs, boolean secondaryReads) {
            super(header, writer, genomeLocParser, referenceReader, maxRecordsInMemory, maxMismatchesInOverhangs, maxBasesInOverhangs, doNotFixOverhangs, secondaryReads);
        }
        @Override
        void fixSplit(SplitRead read, Splice splice) {
            if ( read.read.getStart() == 10000) {
                final GATKRead clippedRead = ReadClipper.softClipByReadCoordinates(read.read, 0, 2);
                read.setRead(clippedRead);
            }
        }
    }

    @Test
    public void testUnalignedReadNotClearingReads() {
        final OverhangFixingManager manager = new OverhangFixingManager(getHG19Header(), null, hg19GenomeLocParser, hg19ReferenceReader, 100, 1, 30, false, true);
        manager.addSplicePosition("1",2,3);
        Assert.assertEquals(manager.getReadsInQueueForTesting().size(), 0);
        manager.addReadGroup(Collections.singletonList(ArtificialReadUtils.createArtificialRead(hg19Header, "read1", 1, 10000, new byte[]{(byte)'A'}, new byte[]{(byte)'A'}, "6M")));
        Assert.assertEquals(manager.getReadsInQueueForTesting().size(), 1);
        manager.addReadGroup(Collections.singletonList(ArtificialReadUtils.createArtificialUnmappedRead(hg19Header,new byte[]{(byte)'A'}, new byte[]{(byte)'A'})));
        Assert.assertEquals(manager.getReadsInQueueForTesting().size(), 2);
        manager.addReadGroup(Collections.singletonList(ArtificialReadUtils.createArtificialRead(hg19Header, "read1", 1, 10000, new byte[]{(byte)'A'}, new byte[]{(byte)'A'}, "6M")));
        Assert.assertEquals(manager.getReadsInQueueForTesting().size(), 3);
        // testing that it does clear properly for a new contig
        manager.addReadGroup(Collections.singletonList(ArtificialReadUtils.createArtificialRead(hg19Header, "read1", 2, 10000, new byte[]{(byte)'A'}, new byte[]{(byte)'A'}, "6M")));
        Assert.assertEquals(manager.getReadsInQueueForTesting().size(), 1);
    }

    @Test
    public void testReadWithDeletionsWindowSpanning() {
        final OverhangFixingManager manager = new OverhangFixingManager(getHG19Header(), null, hg19GenomeLocParser, hg19ReferenceReader, 100, 100, 30, false, true);
        // Create a splice that is going to overlap into the deletion of our read, forcing us to check for mismatches to the reference
        OverhangFixingManager.Splice splice = manager.addSplicePosition("1",6816, 11247);
        // Create a read that has a long deletion relative to its remaining mapped bases after splitting
        GATKRead read = ArtificialReadUtils.createArtificialRead(hg19Header, "read1", 0, 11244, new byte[100], new byte[100], "97S2M18D1M");
        OverhangFixingManager.SplitRead split = manager.getSplitRead(read);
        manager.fixSplit(split,splice);
        // Assert that no splitting happened (and no array exception) by asserting a copy of the read was not placed in split.read
        Assert.assertTrue(split.read==read);
    }

    @Test
    public void testReadWithOnlyInsertionInsideSpan() {
        final OverhangFixingManager manager = new OverhangFixingManager(getHG19Header(), null, hg19GenomeLocParser, hg19ReferenceReader, 100, 100, 30, false, true);
        // Create a splice that is going to overlap into the deletion of our read, forcing us to check for mismatches to the reference
        OverhangFixingManager.Splice splice = manager.addSplicePosition("1",6816, 11247);
        // Create a read that is entirely an insertion inside of the splice to demonstrate it is handled without creating an invalid loc
        GATKRead read = ArtificialReadUtils.createArtificialRead(hg19Header, "read1", 0, 11244, new byte[100], new byte[100], "100I");
        OverhangFixingManager.SplitRead split = manager.getSplitRead(read);
        manager.fixSplit(split,splice);
        // Assert that no splitting happened (and no array exception) by asserting a copy of the read was not placed in split.read
        Assert.assertTrue(split.read==read);
    }


}
