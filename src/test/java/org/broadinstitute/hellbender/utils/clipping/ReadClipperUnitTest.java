package org.broadinstitute.hellbender.utils.clipping;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.ReadClipperTestUtils;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.broadinstitute.hellbender.utils.read.ReadUtils.getSoftEnd;
import static org.broadinstitute.hellbender.utils.read.ReadUtils.getSoftStart;

public final class ReadClipperUnitTest extends GATKBaseTest {
    private List<Cigar> cigarList;
    // 9 is the minimum necessary number to try all combinations of cigar types with guarantee of clipping an element with length = 2
    // and there are already 3 cigar elements on the basic cigar
    private final int maximumCigarElements = 6;

    @BeforeClass
    public void init() {
        cigarList = ReadClipperTestUtils.generateCigarList(maximumCigarElements);
        final Cigar additionalCigar = TextCigarCodec.decode("2M3I5M");
        cigarList.add(additionalCigar);
    }

    @Test
    public void testHardClipBothEndsByReferenceCoordinates() {
        for (final Cigar cigar : cigarList) {
            final GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            final int alnStart = read.getStart();
            final int alnEnd = read.getEnd();
            final int readLength = alnStart - alnEnd;
            for (int i = 0; i < readLength / 2; i++) {
                final GATKRead clippedRead = ReadClipper.hardClipBothEndsByReferenceCoordinates(read, alnStart + i, alnEnd - i);
                Assert.assertTrue(clippedRead.getStart() >= alnStart + i, String.format("Clipped alignment start is less than original read (minus %d): %s -> %s", i, read.getCigar().toString(), clippedRead.getCigar().toString()));
                Assert.assertTrue(clippedRead.getEnd() <= alnEnd + i, String.format("Clipped alignment end is greater than original read (minus %d): %s -> %s", i, read.getCigar().toString(), clippedRead.getCigar().toString()));
                assertUnclippedLimits(read, clippedRead);
            }
        }
    }

    @Test
    public void testHardClipByReferenceCoordinates() {
        for (final Cigar cigar : cigarList) {
            final GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            final int start = getSoftStart(read);
            final int stop = getSoftEnd(read);

            for (int i = start; i <= stop; i++) {
                final GATKRead clipLeft = (new ReadClipper(read)).hardClipByReferenceCoordinates(-1, i);
                if (!clipLeft.isEmpty()) {
                    Assert.assertTrue(clipLeft.getStart() >= Math.min(read.getEnd(), i + 1), String.format("Clipped alignment start (%d) is less the expected (%d): %s -> %s", clipLeft.getStart(), i + 1, read.getCigar().toString(), clipLeft.getCigar().toString()));
                    assertRefAlignmentConsistent(clipLeft);
                    assertReadLengthConsistent(clipLeft);
                }

                final GATKRead clipRight = (new ReadClipper(read)).hardClipByReferenceCoordinates(i, -1);
                if (!clipRight.isEmpty() && clipRight.getStart() <= clipRight.getEnd()) {             // alnStart > alnEnd if the entire read is a soft clip now. We can't test those.
                    Assert.assertTrue(clipRight.getEnd() <= Math.max(read.getStart(), i - 1), String.format("Clipped alignment end (%d) is greater than expected (%d): %s -> %s", clipRight.getEnd(), i - 1, read.getCigar().toString(), clipRight.getCigar().toString()));
                    assertRefAlignmentConsistent(clipRight);
                    assertReadLengthConsistent(clipRight);
                }
            }
        }
    }

    @Test
    public void testHardClipByReferenceCoordinatesLeftTail() {
        for (final Cigar cigar : cigarList) {
            final GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            final int alnStart = read.getStart();
            final int alnEnd = read.getEnd();
            if (getSoftStart(read) == alnStart) {                                                                      // we can't test left clipping if the read has hanging soft clips on the left side
                for (int i = alnStart; i <= alnEnd; i++) {
                    final GATKRead clipLeft = ReadClipper.hardClipByReferenceCoordinatesLeftTail(read, i);

                    if (!clipLeft.isEmpty()) {
                        Assert.assertTrue(clipLeft.getStart() >= i + 1, String.format("Clipped alignment start (%d) is less the expected (%d): %s -> %s", clipLeft.getStart(), i + 1, read.getCigar().toString(), clipLeft.getCigar().toString()));
                        assertRefAlignmentConsistent(clipLeft);
                        assertReadLengthConsistent(clipLeft);
                    }
                }
            }
        }
    }

    @Test
    public void testHardClipByReferenceCoordinatesRightTail() {
        for (final Cigar cigar : cigarList) {
            final GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            final int alnStart = read.getStart();
            final int alnEnd = read.getEnd();
            if (getSoftEnd(read) == alnEnd) {                                                                          // we can't test right clipping if the read has hanging soft clips on the right side
                for (int i = alnStart; i <= alnEnd; i++) {
                    final GATKRead clipRight = ReadClipper.hardClipByReferenceCoordinatesRightTail(read, i);
                    if (!clipRight.isEmpty() && clipRight.getStart() <= clipRight.getEnd()) {         // alnStart > alnEnd if the entire read is a soft clip now. We can't test those.
                        Assert.assertTrue(clipRight.getEnd() <= i - 1, String.format("Clipped alignment end (%d) is greater than expected (%d): %s -> %s", clipRight.getEnd(), i - 1, read.getCigar().toString(), clipRight.getCigar().toString()));
                        assertRefAlignmentConsistent(clipRight);
                        assertReadLengthConsistent(clipRight);
                    }
                }
            }
        }
    }

    @Test
    public void testHardClipLowQualEnds() {
        final byte LOW_QUAL = 2;
        final byte HIGH_QUAL = 30;

        /* create a read for every cigar permutation */
        for (final Cigar cigar : cigarList) {
            final GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            final int readLength = read.getLength();
            final byte[] quals = new byte[readLength];

            for (int nLowQualBases = 0; nLowQualBases < readLength; nLowQualBases++) {

                /*  create a read with nLowQualBases in the left tail */
                Arrays.fill(quals, HIGH_QUAL);
                for (int addLeft = 0; addLeft < nLowQualBases; addLeft++)
                    quals[addLeft] = LOW_QUAL;
                read.setBaseQualities(quals);
                final GATKRead clipLeft = ReadClipper.hardClipLowQualEnds(read, LOW_QUAL);
                checkClippedReadsForLowQualEnds(read, clipLeft, LOW_QUAL, nLowQualBases);

                /* create a read with nLowQualBases in the right tail */
                Arrays.fill(quals, HIGH_QUAL);
                for (int addRight = 0; addRight < nLowQualBases; addRight++)
                    quals[readLength - addRight - 1] = LOW_QUAL;
                read.setBaseQualities(quals);
                final GATKRead clipRight = ReadClipper.hardClipLowQualEnds(read, LOW_QUAL);
                checkClippedReadsForLowQualEnds(read, clipRight, LOW_QUAL, nLowQualBases);

                /* create a read with nLowQualBases on both tails */
                if (nLowQualBases <= readLength / 2) {
                    Arrays.fill(quals, HIGH_QUAL);
                    for (int addBoth = 0; addBoth < nLowQualBases; addBoth++) {
                        quals[addBoth] = LOW_QUAL;
                        quals[readLength - addBoth - 1] = LOW_QUAL;
                    }
                    read.setBaseQualities(quals);
                    final GATKRead clipBoth = ReadClipper.hardClipLowQualEnds(read, LOW_QUAL);
                    checkClippedReadsForLowQualEnds(read, clipBoth, LOW_QUAL, 2*nLowQualBases);
                }
            }
        }
    }

    @Test
    public void testHardClipSoftClippedBases() {
        for (final Cigar cigar : cigarList) {
            final GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            final GATKRead clippedRead = ReadClipper.hardClipSoftClippedBases(read);
            final CigarCounter original = new CigarCounter(read);
            final CigarCounter clipped = new CigarCounter(clippedRead);

            assertUnclippedLimits(read, clippedRead);  // Make sure limits haven't changed
            original.assertHardClippingSoftClips(clipped); // Make sure we have only clipped SOFT_CLIPS
        }
    }

    @Test(enabled = false)
    public void testHardClipLeadingInsertions() {
        for (final Cigar cigar : cigarList) {
            if (startsWithInsertion(cigar)) {
                final GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
                final GATKRead clippedRead = ReadClipper.hardClipLeadingInsertions(read);

                assertUnclippedLimits(read, clippedRead);        // Make sure limits haven't changed

                int expectedLength = read.getLength() - leadingCigarElementLength(read.getCigar(), CigarOperator.INSERTION);
                if (cigarHasElementsDifferentThanInsertionsAndHardClips(read.getCigar())) {
                    expectedLength -= leadingCigarElementLength(CigarUtils.invertCigar(read.getCigar()), CigarOperator.INSERTION);
                }
                if (!clippedRead.isEmpty()) {
                    Assert.assertEquals(expectedLength, clippedRead.getLength(), String.format("%s -> %s", read.getCigar().toString(), clippedRead.getCigar().toString()));  // check that everything else is still there
                    Assert.assertFalse(startsWithInsertion(clippedRead.getCigar()));                                                                                   // check that the insertions are gone
                } else
                    Assert.assertEquals(expectedLength, 0, String.format("expected length: %d", expectedLength));                                                      // check that the read was expected to be fully clipped
            }
        }
    }

    @Test
    public void testRevertSoftClippedBases() {
        for (final Cigar cigar : cigarList) {
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
        for (final Cigar cigar : cigarList) {
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
        final List<Object[]> tests = new ArrayList<>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        for ( final int softStart : Arrays.asList(-10, -1, 0) ) {
            for ( final int alignmentStart : Arrays.asList(1, 10) ) {
                tests.add(new Object[]{softStart, alignmentStart});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "RevertSoftClipsBeforeContig")
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

    private void assertNoLowQualBases(final GATKRead read, final byte low_qual) {
        if (!read.isEmpty()) {
            final byte[] quals = read.getBaseQualities();
            for (int i = 0; i < quals.length; i++)
                Assert.assertFalse(quals[i] <= low_qual, String.format("Found low qual (%d) base after hard clipping. Position: %d -- %s", low_qual, i, read.getCigar().toString()));
        }
    }

    private void checkClippedReadsForLowQualEnds(final GATKRead read, final GATKRead clippedRead, final byte lowQual, final int nLowQualBases) {
        assertNoLowQualBases(clippedRead, lowQual);  // Make sure the low qualities are gone
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

    /**
    * Asserts that the length of alignment on the reference is consistent with the CIGAR
    * after clipping
    *
    * @param clippedRead input read
    * */
    private void assertRefAlignmentConsistent(final GATKRead clippedRead) {
        final int cigarRefLength = clippedRead.getCigar().getReferenceLength();
        final int readRefLength = clippedRead.isUnmapped() ? 0 : clippedRead.getLengthOnReference();

        Assert.assertEquals(cigarRefLength, readRefLength);
    }

    /**
     * Asserts that the length of the read is consistent with the CIGAR after clipping
     *
     * @param clippedRead input read
     * */
    private void assertReadLengthConsistent(final GATKRead clippedRead){
        final int cigarReadLength = clippedRead.getCigar().getReadLength();
        final int readReadLength = clippedRead.getLength();
        Assert.assertEquals(cigarReadLength, readReadLength);
    }

    /**
     * Asserts that number of clipped bases from the read is consistent with the requested clipping
     *
     * @param original original read
     * @param clipped clipped read
     * @param clipping number of clipped bases requested
     * */
    private void assertReadClippingConsistent(final GATKRead original,
                                              final GATKRead clipped,
                                              final int clipping) {
        final int clip_diff = original.getLength() - clipped.getLength();
        Assert.assertEquals(clip_diff, clipping);
    }

    private boolean startsWithInsertion(final Cigar cigar) {
        return leadingCigarElementLength(cigar, CigarOperator.INSERTION) > 0;
    }

    private int leadingCigarElementLength(final Cigar cigar, final CigarOperator operator) {
        for (final CigarElement cigarElement : cigar.getCigarElements()) {
            if (cigarElement.getOperator() == operator)
                return cigarElement.getLength();
            if (cigarElement.getOperator() != CigarOperator.HARD_CLIP)
                break;
        }
        return 0;
    }

    private boolean cigarHasElementsDifferentThanInsertionsAndHardClips(final Cigar cigar) {
        for (final CigarElement cigarElement : cigar.getCigarElements())
            if (cigarElement.getOperator() != CigarOperator.INSERTION && cigarElement.getOperator() != CigarOperator.HARD_CLIP)
                return true;
        return false;
    }

    private class CigarCounter {
        final private Map<CigarOperator, Integer> counter;

        Integer getCounterForOp(final CigarOperator operator) {
            return counter.get(operator);
        }

        CigarCounter(final GATKRead read) {
            final CigarOperator[] operators = CigarOperator.values();
            counter = new LinkedHashMap<>(operators.length);

            for (final CigarOperator op : operators)
                counter.put(op, 0);

            for (final CigarElement cigarElement : read.getCigar().getCigarElements())
                counter.put(cigarElement.getOperator(), counter.get(cigarElement.getOperator()) + cigarElement.getLength());
        }

        boolean assertHardClippingSoftClips(final CigarCounter clipped) {
            for (final CigarOperator op : counter.keySet()) {
                if (op == CigarOperator.HARD_CLIP || op == CigarOperator.SOFT_CLIP) {
                    final int counterTotal = counter.get(CigarOperator.HARD_CLIP) + counter.get(CigarOperator.SOFT_CLIP);
                    final int clippedHard = clipped.getCounterForOp(CigarOperator.HARD_CLIP);
                    final int clippedSoft = clipped.getCounterForOp(CigarOperator.SOFT_CLIP);

                    Assert.assertEquals(counterTotal, clippedHard);
                    Assert.assertEquals(clippedSoft, 0);
                } else {
                    Assert.assertEquals(counter.get(op), clipped.getCounterForOp(op));
                }
            }
            return true;
        }

    }

    @Test
    public void testRevertEntirelySoftClippedReads() {
        final GATKRead read = ReadClipperTestUtils.makeReadFromCigar("2H1S3H");
        final GATKRead clippedRead = ReadClipper.revertSoftClippedBases(read);
        Assert.assertEquals(clippedRead.getStart(), getSoftStart(read));
    }


    @Test
    public void testSoftClipBothEndsByReferenceCoordinates() {
        for (final Cigar cigar : cigarList) {
            final GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            final int alnStart = read.getStart();
            final int alnEnd = read.getEnd();
            final int readLength = alnStart - alnEnd;
            for (int i = 0; i < readLength / 2; i++) {
                final GATKRead clippedRead = ReadClipper.softClipBothEndsByReferenceCoordinates(read, alnStart + i, alnEnd - i);
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

        final GATKRead read = ReadClipperTestUtils.makeReadFromCigar("8M");
        final ReadClipper clipper = new ReadClipper(read);
        final ClippingOp op = new ClippingOp(0, 7);
        clipper.addOp(op);
        final GATKRead softResult = clipper.clipRead(ClippingRepresentation.SOFTCLIP_BASES);
        Assert.assertEquals(softResult.getCigar().toString(), "7S1M");
    }



    //Test pending resolution of issue #2022
    @Test (enabled = false)
    public void testSoftClipByReferenceCoordinates() {
        for (final Cigar cigar : cigarList) {
            if(cigar.isValid(null, -1) != null) {
                continue;
            }
            final GATKRead read = ReadClipperTestUtils.makeReadFromCigar(cigar);
            final int start = getSoftStart(read);
            final int stop = getSoftEnd(read);

            for (int i = start; i <= stop; i++) {
                final GATKRead clipLeft = (new ReadClipper(read.copy())).softClipByReferenceCoordinates(-1, i);
                if (!clipLeft.isEmpty()) {
                    Assert.assertTrue(clipLeft.getStart() >= Math.min(read.getEnd(), i + 1), String.format("Clipped alignment start (%d) is less the expected (%d): %s -> %s", clipLeft.getStart(), i + 1, read.getCigar().toString(), clipLeft.getCigar().toString()));
                }
                final GATKRead clipRight = (new ReadClipper(read.copy())).softClipByReferenceCoordinates(i, -1);
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
