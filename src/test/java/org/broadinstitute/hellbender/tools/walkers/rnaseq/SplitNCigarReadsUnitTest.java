package org.broadinstitute.hellbender.tools.walkers.rnaseq;

import com.google.common.primitives.Ints;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GATKReadWriter;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.testutils.ReadClipperTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * Tests all possible (and valid) cigar strings that might contain any cigar elements. It uses a code that were written to test the ReadClipper walker.
 * For valid cigar sting in length 8 there are few thousands options, with N in every possible option and with more than one N (for example 1M1N1M1N1M1N2M).
 * The cigarElements array is used to provide all the possible cigar element that might be included.
 */
public final class SplitNCigarReadsUnitTest extends GATKBaseTest {
    final static CigarElement[] cigarElements = {
            new CigarElement(1, CigarOperator.HARD_CLIP),
            new CigarElement(1, CigarOperator.SOFT_CLIP),
            new CigarElement(1, CigarOperator.INSERTION),
            new CigarElement(1, CigarOperator.DELETION),
            new CigarElement(1, CigarOperator.MATCH_OR_MISMATCH),
            new CigarElement(1, CigarOperator.SKIPPED_REGION)
    };
    private SAMFileHeader header = new SAMFileHeader();;

    private final class TestManager extends OverhangFixingManager {
        public TestManager( final SAMFileHeader header , DummyTestWriter writer) {
            super(header, writer, hg19GenomeLocParser, hg19ReferenceReader, 10000, 1, 40, false, true);
        }
    }

    @Test
    public void testEmptyReads() {
        DummyTestWriter writer = new DummyTestWriter();
        header.setSequenceDictionary(hg19GenomeLocParser.getSequenceDictionary());
        TestManager manager = new TestManager(header, writer);
        manager.activateWriting();

        //Testing that bogus splits make it through unaffected
        GATKRead read1 = ReadUtils.emptyRead(ReadClipperTestUtils.makeReadFromCigar("1S4N2S3N4H"));
        SplitNCigarReads.splitNCigarRead(read1, manager, true, header, true);
        manager.flush();
        Assert.assertEquals(1, writer.writtenReads.size());
        Assert.assertEquals("*", writer.writtenReads.get(0).getCigar().toString());
    }

    @Test
    public void testBogusNSplits() {
        DummyTestWriter writer = new DummyTestWriter();
        header.setSequenceDictionary(hg19GenomeLocParser.getSequenceDictionary());
        TestManager manager = new TestManager(header, writer);
        manager.activateWriting();

        //Testing that bogus splits make it through unaffected
        GATKRead read1 = ReadClipperTestUtils.makeReadFromCigar("1S4N2S3N4H");
        SplitNCigarReads.splitNCigarRead(read1, manager, true, header, true);
        manager.flush();
        Assert.assertEquals(1, writer.writtenReads.size());
        Assert.assertEquals("1S4N2S3N4H", writer.writtenReads.get(0).getCigar().toString());
    }

    @Test
    public void testBogusMidNSection() {
        //Testing that bogus subsections dont end up clipped
        DummyTestWriter writer = new DummyTestWriter();
        header.setSequenceDictionary(hg19GenomeLocParser.getSequenceDictionary());
        TestManager manager = new TestManager(header, writer);
        manager.activateWriting();

        manager = new TestManager(header, writer);
        manager.activateWriting();
        GATKRead read2 = ReadClipperTestUtils.makeReadFromCigar("1S3N2M10N1M4H");
        SplitNCigarReads.splitNCigarRead(read2, manager, true, header, true);
        manager.flush();
        Assert.assertEquals(2, writer.writtenReads.size());
        Assert.assertEquals("1S2M1S4H", writer.writtenReads.get(0).getCigar().toString());
        Assert.assertEquals("3S1M4H", writer.writtenReads.get(1).getCigar().toString());
    }

    @Test
    public void testSplitNComplexCase() {
        //Complex Case involving insertions & deletions
        DummyTestWriter writer = new DummyTestWriter();
        header.setSequenceDictionary(hg19GenomeLocParser.getSequenceDictionary());
        TestManager manager = new TestManager(header, writer);

        manager.activateWriting();
        writer = new DummyTestWriter();
        manager = new TestManager(header, writer);
        manager.activateWriting();
        GATKRead read3 = ReadClipperTestUtils.makeReadFromCigar("1H2M2D1M2N1M2I1N1M2S1N2S1N2S");
        read3.setAttribute("MC", "11M4N11M");
        SplitNCigarReads.splitNCigarRead(read3, manager,true, header, true);
        manager.flush();
        Assert.assertEquals(3, writer.writtenReads.size());
        Assert.assertEquals("1H2M2D1M10S", writer.writtenReads.get(0).getCigar().toString());
        Assert.assertEquals("1H3S1M2I7S", writer.writtenReads.get(1).getCigar().toString());
        Assert.assertEquals("1H6S1M6S", writer.writtenReads.get(2).getCigar().toString());
        // Testing that the mate cigar string setting is correct
        Assert.assertEquals("11M11S", writer.writtenReads.get(0).getAttributeAsString("MC"));
        Assert.assertEquals("11M11S", writer.writtenReads.get(1).getAttributeAsString("MC"));
        Assert.assertEquals("11M11S", writer.writtenReads.get(2).getAttributeAsString("MC"));
    }


    @Test
    public void splitReadAtN() {
        final int maxCigarElements = 9;
        final List<Cigar> cigarList = ReadClipperTestUtils.generateCigarList(maxCigarElements, cigarElements);

        // For Debugging use those lines (instead of above cigarList) to create specific read:
        //------------------------------------------------------------------------------------
        // final SAMRecord tmpRead = SAMRecord.createRandomRead(6);
        // tmpRead.setCigarString("1M1N1M");

        // final List<Cigar> cigarList = new ArrayList<>();
        // cigarList.add(tmpRead.getCigar());

        for(Cigar cigar: cigarList){

            final int numOfSplits = numOfNElements(cigar.getCigarElements());

            if(numOfSplits != 0 && isCigarDoesNotHaveEmptyRegionsBetweenNs(cigar)){
                final SAMFileHeader header = new SAMFileHeader();
                header.setSequenceDictionary(hg19GenomeLocParser.getSequenceDictionary());
                final TestManager manager = new TestManager(header, new DummyTestWriter());
                manager.activateWriting();

                GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
                SplitNCigarReads.splitNCigarRead(read, manager, true, header, true);
                List<List<OverhangFixingManager.SplitRead>> splitReadGroups = manager.getReadsInQueueForTesting();
                final int expectedReads = numOfSplits+1;
                Assert.assertEquals(manager.getNReadsInQueue(), expectedReads, "wrong number of reads after split read with cigar: " + cigar + " at Ns [expected]: " + expectedReads + " [actual value]: " + manager.getNReadsInQueue());
                final int readLengths = consecutiveNonNElements(read.getCigar().getCigarElements());
                int index = 0;
                for(final OverhangFixingManager.SplitRead splitRead: splitReadGroups.get(0)){
                    Assert.assertTrue(splitRead.read.getLength() == readLengths,
                            "the " + index + " (starting with 0) split read has a wrong length.\n" +
                                    "cigar of original read: " + cigar + "\n" +
                                    "expected length: " + readLengths + "\n" +
                                    "actual length: " + splitRead.read.getLength() + "\n");
                    assertBases(splitRead.read.getBases(), read.getBases());
                    index++;
                }
            }
        }
    }

    private int numOfNElements(final List<CigarElement> cigarElements){
        return Ints.checkedCast(cigarElements.stream()
                .filter(e -> e.getOperator() == CigarOperator.SKIPPED_REGION)
                .count());
    }

    private static boolean isCigarDoesNotHaveEmptyRegionsBetweenNs(final Cigar cigar) {
        boolean sawM = false;

        for (CigarElement cigarElement : cigar.getCigarElements()) {
            if (cigarElement.getOperator().equals(CigarOperator.SKIPPED_REGION)) {
                if(!sawM)
                    return false;
                sawM = false;
            }
            if (cigarElement.getOperator().equals(CigarOperator.MATCH_OR_MISMATCH))
                sawM = true;

        }
        return sawM;
    }

    private int consecutiveNonNElements(final List<CigarElement> cigarElements){
        int consecutiveLength = 0;
        for(CigarElement element: cigarElements){
            final CigarOperator op = element.getOperator();
            if(op.equals(CigarOperator.MATCH_OR_MISMATCH) || op.equals(CigarOperator.SOFT_CLIP) || op.equals(CigarOperator.INSERTION)){
                consecutiveLength += element.getLength();
            }
        }
        return consecutiveLength;
    }

    private void assertBases(final byte[] actualBase, final byte[] expectedBase) {
        for (int i = 0; i < actualBase.length; i++) {
            Assert.assertEquals(actualBase[i], expectedBase[i], "unmatched bases between: " + Arrays.toString(actualBase) + "\nand:\n" + Arrays.toString(expectedBase) + "\nat position: " + i);
        }
    }

    private static class DummyTestWriter implements GATKReadWriter {
        public List<GATKRead> writtenReads = new ArrayList<>();
        @Override
        public void close() {}
        @Override
        public void addRead(GATKRead read) {
            writtenReads.add(read);
        }
    }
}
