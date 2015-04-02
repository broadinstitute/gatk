package org.broadinstitute.hellbender.tools.walkers.rnaseq;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class OverhangFixingManagerUnitTest extends BaseTest {

    private SAMFileHeader getHG19Header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(hg19GenomeLocParser.getSequenceDictionary());
        return header;
    }

    @Test
    public void testCleanSplices() {

        final OverhangFixingManager manager = new OverhangFixingManager(getHG19Header(), null, hg19GenomeLocParser, hg19ReferenceReader, 10000, 1, 40, false);

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
        final OverhangFixingManager manager = new OverhangFixingManager(getHG19Header(), null, hg19GenomeLocParser, hg19ReferenceReader, 10000, 1, 40, false);
        Assert.assertFalse(manager.overhangingBasesMismatch(read, readStart, ref, refStart, overhang));
    }

    @DataProvider(name = "MismatchTest")
    public Object[][] makeMismatchData() {
        final List<Object[]> tests = new ArrayList<>();

        final byte[] AAAA = new byte[]{(byte)'A', (byte)'A', (byte)'A', (byte)'A'};
        final byte[] AAAC = new byte[]{(byte)'A', (byte)'A', (byte)'A', (byte)'C'};
        final byte[] AAAAAA = new byte[]{(byte)'A', (byte)'A', (byte)'A', (byte)'A', (byte)'A', (byte)'A'};
        final byte[] AAAACA = new byte[]{(byte)'A', (byte)'A', (byte)'A', (byte)'A', (byte)'C', (byte)'A'};
        final byte[] AAAACC = new byte[]{(byte)'A', (byte)'A', (byte)'A', (byte)'A', (byte)'C', (byte)'C'};

        tests.add(new Object[]{AAAA, 2, AAAA, 2, 2, false});
        tests.add(new Object[]{AAAA, 2, AAAC, 2, 2, true});
        tests.add(new Object[]{AAAAAA, 3, AAAACA, 3, 3, false});
        tests.add(new Object[]{AAAAAA, 3, AAAACC, 3, 3, true});
        tests.add(new Object[]{AAAAAA, 4, AAAACC, 4, 2, true});
        tests.add(new Object[]{AAAAAA, 2, AAAACC, 2, 3, false});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MismatchTest")
    public void testMismatch(final byte[] read, final int readStart, final byte[] ref, final int refStart, final int overhang, final boolean expected) {
        final OverhangFixingManager manager = new OverhangFixingManager(getHG19Header(), null, hg19GenomeLocParser, hg19ReferenceReader, 10000, 1, 40, false);
        Assert.assertEquals(manager.overhangingBasesMismatch(read, readStart, ref, refStart, overhang), expected, new String(read) + " vs. " + new String(ref) + " @" + overhang);
    }

    @Test
    public void testUnmappedReadsDoNotFail() {
        // create an unmapped read
        final GATKRead read = ArtificialReadUtils.createRandomRead(100);
        read.setName("foo");
        read.setCigar("*");
        read.setIsUnmapped();

        // try to add it to the manager
        final OverhangFixingManager manager = new OverhangFixingManager(getHG19Header(), null, null, null, 100, 1, 30, false);
        manager.addRead(read); // we just want to make sure that the following call does not fail
        Assert.assertTrue(true);
    }
}
