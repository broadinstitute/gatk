package org.broadinstitute.hellbender.utils.clipping;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import static org.broadinstitute.hellbender.utils.read.ReadUtils.*;

public final class ReadClipperUnitTest extends BaseTest {
    List<Cigar> cigarList;
    int maximumCigarSize = 10;                                                                                           // 6 is the minimum necessary number to try all combinations of cigar types with guarantee of clipping an element with length = 2

    @BeforeClass
    public void init() {
        cigarList = ReadClipperTestUtils.generateCigarList(maximumCigarSize);
    }

    @Test
    public void testHardClipBothEndsByReferenceCoordinates() {
        for (Cigar cigar : cigarList) {
            SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int alnStart = read.getAlignmentStart();
            int alnEnd = read.getAlignmentEnd();
            int readLength = alnStart - alnEnd;
            for (int i = 0; i < readLength / 2; i++) {
                SAMRecord clippedRead = ReadClipper.hardClipBothEndsByReferenceCoordinates(read, alnStart + i, alnEnd - i);
                Assert.assertTrue(clippedRead.getAlignmentStart() >= alnStart + i, String.format("Clipped alignment start is less than original read (minus %d): %s -> %s", i, read.getCigarString(), clippedRead.getCigarString()));
                Assert.assertTrue(clippedRead.getAlignmentEnd() <= alnEnd + i, String.format("Clipped alignment end is greater than original read (minus %d): %s -> %s", i, read.getCigarString(), clippedRead.getCigarString()));
                assertUnclippedLimits(read, clippedRead);
            }
        }
    }

    @Test
    public void testHardClipByReadCoordinates() {
        for (Cigar cigar : cigarList) {
            SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int readLength = read.getReadLength();
            for (int i = 0; i < readLength; i++) {
                SAMRecord clipLeft = ReadClipper.hardClipByReadCoordinates(read, 0, i);
                Assert.assertTrue(clipLeft.getReadLength() <= readLength - i, String.format("Clipped read length is greater than original read length (minus %d): %s -> %s", i, read.getCigarString(), clipLeft.getCigarString()));
                assertUnclippedLimits(read, clipLeft);

                SAMRecord clipRight = ReadClipper.hardClipByReadCoordinates(read, i, readLength - 1);
                Assert.assertTrue(clipRight.getReadLength() <= i, String.format("Clipped read length is greater than original read length (minus %d): %s -> %s", i, read.getCigarString(), clipRight.getCigarString()));
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
        SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(originalReadLength + "M");
        read.getReadLength(); // provoke the caching of the read length
        final int expectedReadLength = originalReadLength - nToClip;
        SAMRecord clipped = ReadClipper.hardClipByReadCoordinates(read, 0, nToClip - 1);
        Assert.assertEquals(clipped.getReadLength(), expectedReadLength,
                String.format("Clipped read length %d with cigar %s not equal to the expected read length %d after clipping %d bases from the left from a %d bp read with cigar %s",
                        clipped.getReadLength(), clipped.getCigar(), expectedReadLength, nToClip, read.getReadLength(), read.getCigar()));
    }

    @Test
    public void testHardClipByReferenceCoordinates() {
        for (Cigar cigar : cigarList) {
            SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int start = getSoftStart(read);
            int stop = getSoftEnd(read);

            for (int i = start; i <= stop; i++) {
                SAMRecord clipLeft = (new ReadClipper(read)).hardClipByReferenceCoordinates(-1, i);
                if (!isEmpty(clipLeft)) {
                    Assert.assertTrue(clipLeft.getAlignmentStart() >= Math.min(read.getAlignmentEnd(), i + 1), String.format("Clipped alignment start (%d) is less the expected (%d): %s -> %s", clipLeft.getAlignmentStart(), i + 1, read.getCigarString(), clipLeft.getCigarString()));
                    assertUnclippedLimits(read, clipLeft);
                }

                SAMRecord clipRight = (new ReadClipper(read)).hardClipByReferenceCoordinates(i, -1);
                if (!isEmpty(clipRight) && clipRight.getAlignmentStart() <= clipRight.getAlignmentEnd()) {             // alnStart > alnEnd if the entire read is a soft clip now. We can't test those.
                    Assert.assertTrue(clipRight.getAlignmentEnd() <= Math.max(read.getAlignmentStart(), i - 1), String.format("Clipped alignment end (%d) is greater than expected (%d): %s -> %s", clipRight.getAlignmentEnd(), i - 1, read.getCigarString(), clipRight.getCigarString()));
                    assertUnclippedLimits(read, clipRight);
                }
            }
        }
    }

    @Test
    public void testHardClipByReferenceCoordinatesLeftTail() {
        for (Cigar cigar : cigarList) {
            SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int alnStart = read.getAlignmentStart();
            int alnEnd = read.getAlignmentEnd();
            if (getSoftStart(read) == alnStart) {                                                                      // we can't test left clipping if the read has hanging soft clips on the left side
                for (int i = alnStart; i <= alnEnd; i++) {
                    SAMRecord clipLeft = ReadClipper.hardClipByReferenceCoordinatesLeftTail(read, i);

                    if (!isEmpty(clipLeft)) {
                        Assert.assertTrue(clipLeft.getAlignmentStart() >= i + 1, String.format("Clipped alignment start (%d) is less the expected (%d): %s -> %s", clipLeft.getAlignmentStart(), i + 1, read.getCigarString(), clipLeft.getCigarString()));
                        assertUnclippedLimits(read, clipLeft);
                    }
                }
            }
        }
    }

    @Test
    public void testHardClipByReferenceCoordinatesRightTail() {
        for (Cigar cigar : cigarList) {
            SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int alnStart = read.getAlignmentStart();
            int alnEnd = read.getAlignmentEnd();
            if (getSoftEnd(read) == alnEnd) {                                                                          // we can't test right clipping if the read has hanging soft clips on the right side
                for (int i = alnStart; i <= alnEnd; i++) {
                    SAMRecord clipRight = ReadClipper.hardClipByReferenceCoordinatesRightTail(read, i);
                    if (!isEmpty(clipRight) && clipRight.getAlignmentStart() <= clipRight.getAlignmentEnd()) {         // alnStart > alnEnd if the entire read is a soft clip now. We can't test those.
                        Assert.assertTrue(clipRight.getAlignmentEnd() <= i - 1, String.format("Clipped alignment end (%d) is greater than expected (%d): %s -> %s", clipRight.getAlignmentEnd(), i - 1, read.getCigarString(), clipRight.getCigarString()));
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
            SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            int readLength = read.getReadLength();
            byte[] quals = new byte[readLength];

            for (int nLowQualBases = 0; nLowQualBases < readLength; nLowQualBases++) {

                /**  create a read with nLowQualBases in the left tail */
                Utils.fillArrayWithByte(quals, HIGH_QUAL);
                for (int addLeft = 0; addLeft < nLowQualBases; addLeft++)
                    quals[addLeft] = LOW_QUAL;
                read.setBaseQualities(quals);
                SAMRecord clipLeft = ReadClipper.hardClipLowQualEnds(read, LOW_QUAL);
                checkClippedReadsForLowQualEnds(read, clipLeft, LOW_QUAL, nLowQualBases);

                /** create a read with nLowQualBases in the right tail */
                Utils.fillArrayWithByte(quals, HIGH_QUAL);
                for (int addRight = 0; addRight < nLowQualBases; addRight++)
                    quals[readLength - addRight - 1] = LOW_QUAL;
                read.setBaseQualities(quals);
                SAMRecord clipRight = ReadClipper.hardClipLowQualEnds(read, LOW_QUAL);
                checkClippedReadsForLowQualEnds(read, clipRight, LOW_QUAL, nLowQualBases);

                /** create a read with nLowQualBases on both tails */
                if (nLowQualBases <= readLength / 2) {
                    Utils.fillArrayWithByte(quals, HIGH_QUAL);
                    for (int addBoth = 0; addBoth < nLowQualBases; addBoth++) {
                        quals[addBoth] = LOW_QUAL;
                        quals[readLength - addBoth - 1] = LOW_QUAL;
                    }
                    read.setBaseQualities(quals);
                    SAMRecord clipBoth = ReadClipper.hardClipLowQualEnds(read, LOW_QUAL);
                    checkClippedReadsForLowQualEnds(read, clipBoth, LOW_QUAL, 2*nLowQualBases);
                }
            }
        }
    }

    @Test
    public void testHardClipSoftClippedBases() {
        for (Cigar cigar : cigarList) {
            SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            SAMRecord clippedRead = ReadClipper.hardClipSoftClippedBases(read);
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
                SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
                SAMRecord clippedRead = ReadClipper.hardClipLeadingInsertions(read);

                assertUnclippedLimits(read, clippedRead);        // Make sure limits haven't changed

                int expectedLength = read.getReadLength() - leadingCigarElementLength(read.getCigar(), CigarOperator.INSERTION);
                if (cigarHasElementsDifferentThanInsertionsAndHardClips(read.getCigar()))
                    expectedLength -= leadingCigarElementLength(CigarUtils.invertCigar(read.getCigar()), CigarOperator.INSERTION);

                if (!isEmpty(clippedRead)) {
                    Assert.assertEquals(expectedLength, clippedRead.getReadLength(), String.format("%s -> %s", read.getCigarString(), clippedRead.getCigarString()));  // check that everything else is still there
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

            final SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            final SAMRecord unclipped = ReadClipper.revertSoftClippedBases(read);

            assertUnclippedLimits(read, unclipped);                                                                     // Make sure limits haven't changed

            if (leadingSoftClips > 0 || tailSoftClips > 0) {
                final int expectedStart = read.getAlignmentStart() - leadingSoftClips;
                final int expectedEnd = read.getAlignmentEnd() + tailSoftClips;

                Assert.assertEquals(unclipped.getAlignmentStart(), expectedStart);
                Assert.assertEquals(unclipped.getAlignmentEnd(), expectedEnd);
            } else
                Assert.assertEquals(read.getCigarString(), unclipped.getCigarString());
        }
    }

    @Test
    public void testRevertSoftClippedBasesWithThreshold() {
        for (Cigar cigar : cigarList) {
            final int leadingSoftClips = leadingCigarElementLength(cigar, CigarOperator.SOFT_CLIP);
            final int tailSoftClips = leadingCigarElementLength(CigarUtils.invertCigar(cigar), CigarOperator.SOFT_CLIP);

            final SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            final SAMRecord unclipped = ReadClipper.revertSoftClippedBases(read);

            assertUnclippedLimits(read, unclipped);                                                                     // Make sure limits haven't changed
            Assert.assertNull(read.getCigar().isValid(null, -1));
            Assert.assertNull(unclipped.getCigar().isValid(null, -1));

            if (!(leadingSoftClips > 0 || tailSoftClips > 0))
                Assert.assertEquals(read.getCigarString(), unclipped.getCigarString());

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
        final SAMRecord read = ReadClipperTestUtils.makeReadFromCigar(cigar);
        read.setAlignmentStart(alignmentStart);

        Assert.assertEquals(getSoftStart(read), softStart);
        Assert.assertEquals(read.getAlignmentStart(), alignmentStart);
        Assert.assertEquals(read.getCigarString(), cigar);

        final SAMRecord reverted = ReadClipper.revertSoftClippedBases(read);

        final int expectedAlignmentStart = 1;
        final String expectedCigar = (1 - softStart) + "H" + read.getAlignmentEnd() + "M";
        Assert.assertEquals(getSoftStart(reverted), expectedAlignmentStart);
        Assert.assertEquals(reverted.getAlignmentStart(), expectedAlignmentStart);
        Assert.assertEquals(reverted.getCigarString(), expectedCigar);
    }

    private void assertNoLowQualBases(SAMRecord read, byte low_qual) {
        if (!isEmpty(read)) {
            byte[] quals = read.getBaseQualities();
            for (int i = 0; i < quals.length; i++)
                Assert.assertFalse(quals[i] <= low_qual, String.format("Found low qual (%d) base after hard clipping. Position: %d -- %s", low_qual, i, read.getCigarString()));
        }
    }

    private void checkClippedReadsForLowQualEnds(SAMRecord read, SAMRecord clippedRead, byte lowQual, int nLowQualBases) {
        assertUnclippedLimits(read, clippedRead);                                                                       // Make sure limits haven't changed
        assertNoLowQualBases(clippedRead, lowQual);                                                                     // Make sure the low qualities are gone
    }

    /**
     * Asserts that clipping doesn't change the getUnclippedStart / getUnclippedEnd
     *
     * @param original original read
     * @param clipped clipped read
     */
    private void assertUnclippedLimits(final SAMRecord original, final SAMRecord clipped) {
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
        private HashMap<CigarOperator, Integer> counter;

        public Integer getCounterForOp(CigarOperator operator) {
            return counter.get(operator);
        }

        public CigarCounter(SAMRecord read) {
            CigarOperator[] operators = CigarOperator.values();
            counter = new HashMap<>(operators.length);

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
    public void testRevertEntirelySoftclippedReads() {
        SAMRecord read = ReadClipperTestUtils.makeReadFromCigar("2H1S3H");
        SAMRecord clippedRead = ReadClipper.revertSoftClippedBases(read);
        Assert.assertEquals(clippedRead.getAlignmentStart(), getSoftStart(read));
    }

}