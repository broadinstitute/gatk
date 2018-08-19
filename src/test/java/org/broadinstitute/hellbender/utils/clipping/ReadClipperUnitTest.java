package org.broadinstitute.hellbender.utils.clipping;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.ReadClipperTestUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

import static org.broadinstitute.hellbender.utils.read.ReadUtils.getSoftEnd;
import static org.broadinstitute.hellbender.utils.read.ReadUtils.getSoftStart;

public final class ReadClipperUnitTest extends GATKBaseTest {
    List<Cigar> cigarList;
    int maximumCigarElements = 9;                                                                                           // 6 is the minimum necessary number to try all combinations of cigar types with guarantee of clipping an element with length = 2

    @BeforeClass
    public void init() {
        cigarList = ReadClipperTestUtils.generateCigarList(maximumCigarElements);
    }

    @Test
    public void testHardClipBothEndsByReferenceCoordinates() {
        for (Cigar cigar : cigarList) {
            GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int alnStart = read.getStart();
            int alnEnd = read.getEnd();
            int readLength = alnStart - alnEnd;
            for (int i = 0; i < readLength / 2; i++) {
                GATKRead clippedRead = ReadClipper.hardClipBothEndsByReferenceCoordinates(read, alnStart + i, alnEnd - i);
                Assert.assertTrue(clippedRead.getStart() >= alnStart + i, String.format("Clipped alignment start is less than original read (minus %d): %s -> %s", i, read.getCigar().toString(), clippedRead.getCigar().toString()));
                Assert.assertTrue(clippedRead.getEnd() <= alnEnd + i, String.format("Clipped alignment end is greater than original read (minus %d): %s -> %s", i, read.getCigar().toString(), clippedRead.getCigar().toString()));
                assertUnclippedLimits(read, clippedRead);
            }
        }
    }

    @Test
    public void testHardClipByReadCoordinates() {
        for (Cigar cigar : cigarList) {
            GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int readLength = read.getLength();
            for (int i = 0; i < readLength; i++) {
                GATKRead clipLeft = ReadClipper.hardClipByReadCoordinates(read, 0, i);
                Assert.assertTrue(clipLeft.getLength() <= readLength - i, String.format("Clipped read length is greater than original read length (minus %d): %s -> %s", i, read.getCigar().toString(), clipLeft.getCigar().toString()));
                assertUnclippedLimits(read, clipLeft);

                GATKRead clipRight = ReadClipper.hardClipByReadCoordinates(read, i, readLength - 1);
                Assert.assertTrue(clipRight.getLength() <= i, String.format("Clipped read length is greater than original read length (minus %d): %s -> %s", i, read.getCigar().toString(), clipRight.getCigar().toString()));
                assertUnclippedLimits(read, clipRight);
            }
        }
    }

    @DataProvider(name = "ClippedReadLengthData")
    public Object[][] makeClippedReadLengthData() {
        final List<Object[]> tests = new ArrayList<>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        final int originalReadLength = 50;
        for ( int nToClip = 1; nToClip < originalReadLength - 1; nToClip++ ) {
            tests.add(new Object[]{originalReadLength, nToClip});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ClippedReadLengthData")
    public void testHardClipReadLengthIsRight(final int originalReadLength, final int nToClip) {
        GATKRead read = ReadClipperTestUtils.makeReadFromCigar(originalReadLength + "M");
        read.getLength(); // provoke the caching of the read length
        final int expectedReadLength = originalReadLength - nToClip;
        GATKRead clipped = ReadClipper.hardClipByReadCoordinates(read, 0, nToClip - 1);
        Assert.assertEquals(clipped.getLength(), expectedReadLength,
                String.format("Clipped read length %d with cigar %s not equal to the expected read length %d after clipping %d bases from the left from a %d bp read with cigar %s",
                        clipped.getLength(), clipped.getCigar(), expectedReadLength, nToClip, read.getLength(), read.getCigar()));
    }

    @Test
    public void testHardClipByReferenceCoordinates() {
        for (Cigar cigar : cigarList) {
            GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int start = getSoftStart(read);
            int stop = getSoftEnd(read);

            for (int i = start; i <= stop; i++) {
                GATKRead clipLeft = (new ReadClipper(read)).hardClipByReferenceCoordinates(-1, i);
                if (!clipLeft.isEmpty()) {
                    Assert.assertTrue(clipLeft.getStart() >= Math.min(read.getEnd(), i + 1), String.format("Clipped alignment start (%d) is less the expected (%d): %s -> %s", clipLeft.getStart(), i + 1, read.getCigar().toString(), clipLeft.getCigar().toString()));
                    assertUnclippedLimits(read, clipLeft);
                }

                GATKRead clipRight = (new ReadClipper(read)).hardClipByReferenceCoordinates(i, -1);
                if (!clipRight.isEmpty() && clipRight.getStart() <= clipRight.getEnd()) {             // alnStart > alnEnd if the entire read is a soft clip now. We can't test those.
                    Assert.assertTrue(clipRight.getEnd() <= Math.max(read.getStart(), i - 1), String.format("Clipped alignment end (%d) is greater than expected (%d): %s -> %s", clipRight.getEnd(), i - 1, read.getCigar().toString(), clipRight.getCigar().toString()));
                    assertUnclippedLimits(read, clipRight);
                }
            }
        }
    }

    @Test
    public void testHardClipByReferenceCoordinatesLeftTail() {
        for (Cigar cigar : cigarList) {
            GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int alnStart = read.getStart();
            int alnEnd = read.getEnd();
            if (getSoftStart(read) == alnStart) {                                                                      // we can't test left clipping if the read has hanging soft clips on the left side
                for (int i = alnStart; i <= alnEnd; i++) {
                    GATKRead clipLeft = ReadClipper.hardClipByReferenceCoordinatesLeftTail(read, i);

                    if (!clipLeft.isEmpty()) {
                        Assert.assertTrue(clipLeft.getStart() >= i + 1, String.format("Clipped alignment start (%d) is less the expected (%d): %s -> %s", clipLeft.getStart(), i + 1, read.getCigar().toString(), clipLeft.getCigar().toString()));
                        assertUnclippedLimits(read, clipLeft);
                    }
                }
            }
        }
    }

    @Test
    public void testHardClipByReferenceCoordinatesRightTail() {
        for (Cigar cigar : cigarList) {
            GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int alnStart = read.getStart();
            int alnEnd = read.getEnd();
            if (getSoftEnd(read) == alnEnd) {                                                                          // we can't test right clipping if the read has hanging soft clips on the right side
                for (int i = alnStart; i <= alnEnd; i++) {
                    GATKRead clipRight = ReadClipper.hardClipByReferenceCoordinatesRightTail(read, i);
                    if (!clipRight.isEmpty() && clipRight.getStart() <= clipRight.getEnd()) {         // alnStart > alnEnd if the entire read is a soft clip now. We can't test those.
                        Assert.assertTrue(clipRight.getEnd() <= i - 1, String.format("Clipped alignment end (%d) is greater than expected (%d): %s -> %s", clipRight.getEnd(), i - 1, read.getCigar().toString(), clipRight.getCigar().toString()));
                        assertUnclippedLimits(read, clipRight);
                    }
                }
            }
        }
    }

    @Test
    public void testHardClipLowQualEnds() {
        final byte LOW_QUAL = 2;
        final byte HIGH_QUAL = 30;

        /** create a read for every cigar permutation */
        for (Cigar cigar : cigarList) {
            GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int readLength = read.getLength();
            byte[] quals = new byte[readLength];

            for (int nLowQualBases = 0; nLowQualBases < readLength; nLowQualBases++) {

                /**  create a read with nLowQualBases in the left tail */
                Arrays.fill(quals, HIGH_QUAL);
                for (int addLeft = 0; addLeft < nLowQualBases; addLeft++)
                    quals[addLeft] = LOW_QUAL;
                read.setBaseQualities(quals);
                GATKRead clipLeft = ReadClipper.hardClipLowQualEnds(read, LOW_QUAL);
                checkClippedReadsForLowQualEnds(read, clipLeft, LOW_QUAL, nLowQualBases);

                /** create a read with nLowQualBases in the right tail */
                Arrays.fill(quals, HIGH_QUAL);
                for (int addRight = 0; addRight < nLowQualBases; addRight++)
                    quals[readLength - addRight - 1] = LOW_QUAL;
                read.setBaseQualities(quals);
                GATKRead clipRight = ReadClipper.hardClipLowQualEnds(read, LOW_QUAL);
                checkClippedReadsForLowQualEnds(read, clipRight, LOW_QUAL, nLowQualBases);

                /** create a read with nLowQualBases on both tails */
                if (nLowQualBases <= readLength / 2) {
                    Arrays.fill(quals, HIGH_QUAL);
                    for (int addBoth = 0; addBoth < nLowQualBases; addBoth++) {
                        quals[addBoth] = LOW_QUAL;
                        quals[readLength - addBoth - 1] = LOW_QUAL;
                    }
                    read.setBaseQualities(quals);
                    GATKRead clipBoth = ReadClipper.hardClipLowQualEnds(read, LOW_QUAL);
                    checkClippedReadsForLowQualEnds(read, clipBoth, LOW_QUAL, 2*nLowQualBases);
                }
            }
        }
    }

    @Test
    public void testHardClipSoftClippedBases() {
        for (Cigar cigar : cigarList) {
            GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            GATKRead clippedRead = ReadClipper.hardClipSoftClippedBases(read);
            CigarCounter original = new CigarCounter(read);
            CigarCounter clipped = new CigarCounter(clippedRead);

            assertUnclippedLimits(read, clippedRead);                                                                   // Make sure limits haven't changed
            original.assertHardClippingSoftClips(clipped);                                                              // Make sure we have only clipped SOFT_CLIPS
        }
    }

    @Test(enabled = false)
    public void testHardClipLeadingInsertions() {
        for (Cigar cigar : cigarList) {
            if (startsWithInsertion(cigar)) {
                GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
                GATKRead clippedRead = ReadClipper.hardClipLeadingInsertions(read);

                assertUnclippedLimits(read, clippedRead);        // Make sure limits haven't changed

                int expectedLength = read.getLength() - leadingCigarElementLength(read.getCigar(), CigarOperator.INSERTION);
                if (cigarHasElementsDifferentThanInsertionsAndHardClips(read.getCigar()))
                    expectedLength -= leadingCigarElementLength(CigarUtils.invertCigar(read.getCigar()), CigarOperator.INSERTION);

                if (!clippedRead.isEmpty()) {
                    Assert.assertEquals(expectedLength, clippedRead.getLength(), String.format("%s -> %s", read.getCigar().toString(), clippedRead.getCigar().toString()));  // check that everything else is still there
                    Assert.assertFalse(startsWithInsertion(clippedRead.getCigar()));                                                                                   // check that the insertions are gone
                } else
                    Assert.assertTrue(expectedLength == 0, String.format("expected length: %d", expectedLength));                                                      // check that the read was expected to be fully clipped
            }
        }
    }

    @Test
    public void testRevertSoftClippedBases() {
        for (Cigar cigar : cigarList) {
            final int leadingSoftClips = leadingCigarElementLength(cigar, CigarOperator.SOFT_CLIP);
            final int tailSoftClips = leadingCigarElementLength(CigarUtils.invertCigar(cigar), CigarOperator.SOFT_CLIP);

            final GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            final GATKRead unclipped = ReadClipper.revertSoftClippedBases(read);

            assertUnclippedLimits(read, unclipped);                                                                     // Make sure limits haven't changed

            if (leadingSoftClips > 0 || tailSoftClips > 0) {
                final int expectedStart = read.getStart() - leadingSoftClips;
                final int expectedEnd = read.getEnd() + tailSoftClips;

                Assert.assertEquals(unclipped.getStart(), expectedStart);
                Assert.assertEquals(unclipped.getEnd(), expectedEnd);
            } else
                Assert.assertEquals(read.getCigar().toString(), unclipped.getCigar().toString());
        }
    }

    @Test
    public void testRevertSoftClippedBasesWithThreshold() {
        for (Cigar cigar : cigarList) {
            final int leadingSoftClips = leadingCigarElementLength(cigar, CigarOperator.SOFT_CLIP);
            final int tailSoftClips = leadingCigarElementLength(CigarUtils.invertCigar(cigar), CigarOperator.SOFT_CLIP);

            final GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            final GATKRead unclipped = ReadClipper.revertSoftClippedBases(read);

            assertUnclippedLimits(read, unclipped);                                                                     // Make sure limits haven't changed
            Assert.assertNull(read.getCigar().isValid(null, -1));
            Assert.assertNull(unclipped.getCigar().isValid(null, -1));

            if (!(leadingSoftClips > 0 || tailSoftClips > 0))
                Assert.assertEquals(read.getCigar().toString(), unclipped.getCigar().toString());

        }
    }

    @DataProvider(name = "RevertSoftClipsBeforeContig")
    public Object[][] makeRevertSoftClipsBeforeContig() {
        List<Object[]> tests = new ArrayList<>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        for ( int softStart : Arrays.asList(-10, -1, 0) ) {
            for ( int alignmentStart : Arrays.asList(1, 10) ) {
                tests.add(new Object[]{softStart, alignmentStart});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "RevertSoftClipsBeforeContig")
    public void testRevertSoftClippedBasesBeforeStartOfContig(final int softStart, final int alignmentStart) {
        final int nMatches = 10;
        final int nSoft = -1 * (softStart - alignmentStart);
        final String cigar = nSoft + "S" + nMatches + "M";
        final GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
        read.setPosition(read.getContig(), alignmentStart);

        Assert.assertEquals(getSoftStart(read), softStart);
        Assert.assertEquals(read.getStart(), alignmentStart);
        Assert.assertEquals(read.getCigar().toString(), cigar);

        final GATKRead reverted = ReadClipper.revertSoftClippedBases(read);

        final int expectedAlignmentStart = 1;
        final String expectedCigar = (1 - softStart) + "H" + read.getEnd() + "M";
        Assert.assertEquals(getSoftStart(reverted), expectedAlignmentStart);
        Assert.assertEquals(reverted.getStart(), expectedAlignmentStart);
        Assert.assertEquals(reverted.getCigar().toString(), expectedCigar);
    }

    private void assertNoLowQualBases(GATKRead read, byte low_qual) {
        if (!read.isEmpty()) {
            byte[] quals = read.getBaseQualities();
            for (int i = 0; i < quals.length; i++)
                Assert.assertFalse(quals[i] <= low_qual, String.format("Found low qual (%d) base after hard clipping. Position: %d -- %s", low_qual, i, read.getCigar().toString()));
        }
    }

    private void checkClippedReadsForLowQualEnds(GATKRead read, GATKRead clippedRead, byte lowQual, int nLowQualBases) {
        assertUnclippedLimits(read, clippedRead);                                                                       // Make sure limits haven't changed
        assertNoLowQualBases(clippedRead, lowQual);                                                                     // Make sure the low qualities are gone
    }

    /**
     * Asserts that clipping doesn't change the getUnclippedStart / getUnclippedEnd
     *
     * @param original original read
     * @param clipped clipped read
     */
    private void assertUnclippedLimits(final GATKRead original, final GATKRead clipped) {
        if (CigarUtils.hasNonClippedBases(clipped.getCigar())) {
            Assert.assertEquals(original.getUnclippedStart(), clipped.getUnclippedStart());
            Assert.assertEquals(original.getUnclippedEnd(), clipped.getUnclippedEnd());
        }
    }

    private boolean startsWithInsertion(Cigar cigar) {
        return leadingCigarElementLength(cigar, CigarOperator.INSERTION) > 0;
    }

    private int leadingCigarElementLength(Cigar cigar, CigarOperator operator) {
        for (CigarElement cigarElement : cigar.getCigarElements()) {
            if (cigarElement.getOperator() == operator)
                return cigarElement.getLength();
            if (cigarElement.getOperator() != CigarOperator.HARD_CLIP)
                break;
        }
        return 0;
    }

    private boolean cigarHasElementsDifferentThanInsertionsAndHardClips(Cigar cigar) {
        for (CigarElement cigarElement : cigar.getCigarElements())
            if (cigarElement.getOperator() != CigarOperator.INSERTION && cigarElement.getOperator() != CigarOperator.HARD_CLIP)
                return true;
        return false;
    }

    private class CigarCounter {
        private Map<CigarOperator, Integer> counter;

        public Integer getCounterForOp(CigarOperator operator) {
            return counter.get(operator);
        }

        public CigarCounter(GATKRead read) {
            CigarOperator[] operators = CigarOperator.values();
            counter = new LinkedHashMap<>(operators.length);

            for (CigarOperator op : operators)
                counter.put(op, 0);

            for (CigarElement cigarElement : read.getCigar().getCigarElements())
                counter.put(cigarElement.getOperator(), counter.get(cigarElement.getOperator()) + cigarElement.getLength());
        }

        public boolean assertHardClippingSoftClips(CigarCounter clipped) {
            for (CigarOperator op : counter.keySet()) {
                if (op == CigarOperator.HARD_CLIP || op == CigarOperator.SOFT_CLIP) {
                    int counterTotal = counter.get(CigarOperator.HARD_CLIP) + counter.get(CigarOperator.SOFT_CLIP);
                    int clippedHard = clipped.getCounterForOp(CigarOperator.HARD_CLIP);
                    int clippedSoft = clipped.getCounterForOp(CigarOperator.SOFT_CLIP);

                    Assert.assertEquals(counterTotal, clippedHard);
                    Assert.assertTrue(clippedSoft == 0);
                } else
                    Assert.assertEquals(counter.get(op), clipped.getCounterForOp(op));
            }
            return true;
        }

    }

    @Test
    public void testRevertEntirelySoftClippedReads() {
        GATKRead read = ReadClipperTestUtils.makeReadFromCigar("2H1S3H");
        GATKRead clippedRead = ReadClipper.revertSoftClippedBases(read);
        Assert.assertEquals(clippedRead.getStart(), getSoftStart(read));
    }


    @Test
    public void testSoftClipBothEndsByReferenceCoordinates() {
        for (Cigar cigar : cigarList) {
            GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int alnStart = read.getStart();
            int alnEnd = read.getEnd();
            int readLength = alnStart - alnEnd;
            for (int i = 0; i < readLength / 2; i++) {
                GATKRead clippedRead = ReadClipper.softClipBothEndsByReferenceCoordinates(read, alnStart + i, alnEnd - i);
                Assert.assertTrue(clippedRead.getStart() >= alnStart + i, String.format("Clipped alignment start is less than original read (minus %d): %s -> %s", i, read.getCigar().toString(), clippedRead.getCigar().toString()));
                Assert.assertTrue(clippedRead.getEnd() <= alnEnd + i, String.format("Clipped alignment end is greater than original read (minus %d): %s -> %s", i, read.getCigar().toString(), clippedRead.getCigar().toString()));
                assertUnclippedLimits(read, clippedRead);
            }
        }
    }

    // This test depends on issue #2022 as it tests the current behavior of the clipping operation
    @Test
    public void testSoftClippingOpEdgeCase() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(hg19GenomeLocParser.getSequenceDictionary());

        GATKRead read = ReadClipperTestUtils.makeReadFromCigar("8M");
        ReadClipper clipper = new ReadClipper(read);
        ClippingOp op = new ClippingOp(0, 7);
        clipper.addOp(op);
        GATKRead softResult = clipper.clipRead(ClippingRepresentation.SOFTCLIP_BASES);
        Assert.assertEquals(softResult.getCigar().toString(), "7S1M");
    }



    //Test pending resolution of issue #2022
    @Test (enabled = false)
    public void testSoftClipByReferenceCoordinates() {
        for (Cigar cigar : cigarList) {
            if(cigar.isValid(null, -1) != null) {
                continue;
            }
            GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int start = getSoftStart(read);
            int stop = getSoftEnd(read);

            for (int i = start; i <= stop; i++) {
                GATKRead clipLeft = (new ReadClipper(read.copy())).softClipByReferenceCoordinates(-1, i);
                if (!clipLeft.isEmpty()) {
                    Assert.assertTrue(clipLeft.getStart() >= Math.min(read.getEnd(), i + 1), String.format("Clipped alignment start (%d) is less the expected (%d): %s -> %s", clipLeft.getStart(), i + 1, read.getCigar().toString(), clipLeft.getCigar().toString()));
                }
                GATKRead clipRight = (new ReadClipper(read.copy())).softClipByReferenceCoordinates(i, -1);
                if (!clipRight.isEmpty() && clipRight.getStart() <= clipRight.getEnd()) {             // alnStart > alnEnd if the entire read is a soft clip now. We can't test those.
                    Assert.assertTrue(clipRight.getEnd() <= Math.max(read.getStart(), i - 1), String.format("Clipped alignment end (%d) is greater than expected (%d): %s -> %s", clipRight.getEnd(), i - 1, read.getCigar().toString(), clipRight.getCigar().toString()));
                }
            }
        }
    }

    // Test fix for https://github.com/broadinstitute/gatk/issues/3466
    @Test
    public void testHardClipSoftClippedBasesResultsInEmptyReadDontSetNegativeStartPosition() {
        final GATKRead originalRead = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("170H70S"));
        // It's important that the read be near the start of the contig for this test, to test
        // that we don't attempt to set the read's start position to a negative value during clipping.
        // See https://github.com/broadinstitute/gatk/issues/3466
        originalRead.setPosition(originalRead.getContig(), 100);

        final GATKRead clippedRead = ReadClipper.hardClipSoftClippedBases(originalRead);
        Assert.assertEquals(clippedRead.getLength(), 0);
        Assert.assertTrue(clippedRead.isEmpty());
        Assert.assertEquals(clippedRead.getBases().length, 0);
        Assert.assertEquals(clippedRead.getBaseQualities().length, 0);
        Assert.assertEquals(clippedRead.numCigarElements(), 0);
        Assert.assertTrue(clippedRead.isUnmapped());
    }

    // Test fix for https://github.com/broadinstitute/gatk/issues/3845
    @Test
    public void testRevertSoftClippedBasesDoesntExplodeOnCompletelyClippedRead() {
        final GATKRead originalRead = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("41S59H"));
        // It's important that the read be AT the start of the contig for this test, so that
        // we clip away ALL of the reverted soft-clipped bases, resulting in an empty read.
        originalRead.setPosition(originalRead.getContig(), 1);
        
        final GATKRead clippedRead = ReadClipper.revertSoftClippedBases(originalRead);

        Assert.assertEquals(clippedRead.getLength(), 0);
        Assert.assertTrue(clippedRead.isEmpty());
        Assert.assertEquals(clippedRead.getBases().length, 0);
        Assert.assertEquals(clippedRead.getBaseQualities().length, 0);
        Assert.assertEquals(clippedRead.numCigarElements(), 0);
        Assert.assertTrue(clippedRead.isUnmapped());
    }

}