package org.broadinstitute.hellbender.utils.fragments;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ArtificialSAMUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Test routines for fragments and their collections.
 */
public class FragmentUtilsUnitTest extends BaseTest {
    private static SAMFileHeader header;
    private static SAMReadGroupRecord rgForMerged;

    private class FragmentUtilsTest extends TestDataProvider {
        List<TestState> statesForPileup = new ArrayList<>();
        List<TestState> statesForReads = new ArrayList<>();

        private FragmentUtilsTest(String name, int readLen, int leftStart, int rightStart,
                                  boolean leftIsFirst, boolean leftIsNegative) {
            super(FragmentUtilsTest.class, String.format("%s-leftIsFirst:%b-leftIsNegative:%b", name, leftIsFirst, leftIsNegative));

            List<SAMRecord> pair = ArtificialSAMUtils.createPair(header, "readpair", readLen, leftStart, rightStart, leftIsFirst, leftIsNegative);
            SAMRecord left = pair.get(0);
            SAMRecord right = pair.get(1);

            for ( int pos = leftStart; pos < rightStart + readLen; pos++) {
                boolean posCoveredByLeft = pos >= left.getAlignmentStart() && pos <= left.getAlignmentEnd();
                boolean posCoveredByRight = pos >= right.getAlignmentStart() && pos <= right.getAlignmentEnd();

                if ( posCoveredByLeft || posCoveredByRight ) {
                    List<SAMRecord> reads = new ArrayList<>();
                    List<Integer> offsets = new ArrayList<>();

                    if ( posCoveredByLeft ) {
                        reads.add(left);
                        offsets.add(pos - left.getAlignmentStart());
                    }

                    if ( posCoveredByRight ) {
                        reads.add(right);
                        offsets.add(pos - right.getAlignmentStart());
                    }

                    boolean shouldBeFragment = posCoveredByLeft && posCoveredByRight;
                    ReadPileup pileup = new ReadPileup(null, reads, offsets);
                    TestState testState = new TestState(shouldBeFragment ? 0 : 1, shouldBeFragment ? 1 : 0, pileup, null);
                    statesForPileup.add(testState);
                }

                TestState testState = left.getAlignmentEnd() >= right.getAlignmentStart() ? new TestState(0, 1, null, pair) : new TestState(2, 0, null, pair);
                statesForReads.add(testState);
            }
        }
    }

    private class TestState {
        int expectedSingletons, expectedPairs;
        ReadPileup pileup;
        List<SAMRecord> rawReads;

        private TestState(final int expectedSingletons, final int expectedPairs, final ReadPileup pileup, final List<SAMRecord> rawReads) {
            this.expectedSingletons = expectedSingletons;
            this.expectedPairs = expectedPairs;
            this.pileup = pileup;
            this.rawReads = rawReads;
        }
    }

    @DataProvider(name = "fragmentUtilsTest")
    public Object[][] createTests() {
        for ( boolean leftIsFirst : Arrays.asList(true, false) ) {
            for ( boolean leftIsNegative : Arrays.asList(true, false) ) {
                // Overlapping pair
                // ---->        [first]
                //   <---       [second]
                new FragmentUtilsTest("overlapping-pair", 10, 1, 5, leftIsFirst, leftIsNegative);

                // Non-overlapping pair
                // ---->
                //          <----
                new FragmentUtilsTest("nonoverlapping-pair", 10, 1, 15, leftIsFirst, leftIsNegative);
            }
        }

        return FragmentUtilsTest.getTests(FragmentUtilsTest.class);
    }

    @Test(dataProvider = "fragmentUtilsTest")
    public void testAsPileup(FragmentUtilsTest test) {
        for ( TestState testState : test.statesForPileup ) {
            ReadPileup rbp = testState.pileup;
            FragmentCollection<PileupElement> fp = FragmentCollection.create(rbp);
            Assert.assertEquals(fp.getOverlappingPairs().size(), testState.expectedPairs);
            Assert.assertEquals(fp.getSingletonReads().size(), testState.expectedSingletons);
        }
    }

    @Test(dataProvider = "fragmentUtilsTest")
    public void testAsListOfReadsFromPileup(FragmentUtilsTest test) {
        for ( TestState testState : test.statesForPileup ) {
            FragmentCollection<SAMRecord> fp = FragmentCollection.create(testState.pileup.getReads());
            Assert.assertEquals(fp.getOverlappingPairs().size(), testState.expectedPairs);
            Assert.assertEquals(fp.getSingletonReads().size(), testState.expectedSingletons);
        }
    }

    @Test(dataProvider = "fragmentUtilsTest")
    public void testAsListOfReads(FragmentUtilsTest test) {
        for ( TestState testState : test.statesForReads ) {
            FragmentCollection<SAMRecord> fp = FragmentCollection.create(testState.rawReads);
            Assert.assertEquals(fp.getOverlappingPairs().size(), testState.expectedPairs);
            Assert.assertEquals(fp.getSingletonReads().size(), testState.expectedSingletons);
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testOutOfOrder() {
        final List<SAMRecord> pair = ArtificialSAMUtils.createPair(header, "readpair", 100, 1, 50, true, true);
        final SAMRecord left = pair.get(0);
        final SAMRecord right = pair.get(1);
        final List<SAMRecord> reads = Arrays.asList(right, left); // OUT OF ORDER!
        final List<Integer> offsets = Arrays.asList(0, 50);
        final ReadPileup pileup = new ReadPileup(null, reads, offsets);
        FragmentCollection.create(pileup); // should throw exception
    }

    @BeforeTest
    public void setup() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1,1,1000);
        rgForMerged = new SAMReadGroupRecord("RG1");
    }

    @DataProvider(name = "MergeFragmentsTest")
    public Object[][] createMergeFragmentsTest() throws Exception {
        List<Object[]> tests = new ArrayList<>();

        final String leftFlank = "CCC";
        final String rightFlank = "AAA";
        final String allOverlappingBases = "ACGTACGTGGAACCTTAG";
        for ( int overlapSize = 1; overlapSize < allOverlappingBases.length(); overlapSize++ ) {
            final String overlappingBases = allOverlappingBases.substring(0, overlapSize);
            final byte[] overlappingBaseQuals = new byte[overlapSize];
            for ( int i = 0; i < overlapSize; i++ ) overlappingBaseQuals[i] = (byte)(i + 30);
            final SAMRecord read1  = makeOverlappingRead(leftFlank, 20, overlappingBases, overlappingBaseQuals, "", 30, 1);
            final SAMRecord read2  = makeOverlappingRead("", 20, overlappingBases, overlappingBaseQuals, rightFlank, 30, leftFlank.length() + 1);
            final SAMRecord merged = makeOverlappingRead(leftFlank, 20, overlappingBases, overlappingBaseQuals, rightFlank, 30, 1);
            tests.add(new Object[]{"equalQuals", read1, read2, merged});

            // test that the merged read base quality is the
            tests.add(new Object[]{"lowQualLeft", modifyBaseQualities(read1, leftFlank.length(), overlapSize), read2, merged});
            tests.add(new Object[]{"lowQualRight", read1, modifyBaseQualities(read2, 0, overlapSize), merged});
        }

        return tests.toArray(new Object[][]{});
    }

    private SAMRecord modifyBaseQualities(final SAMRecord read, final int startOffset, final int length) throws Exception {
        final SAMRecord readWithLowQuals = (SAMRecord)read.clone();
        final byte[] withLowQuals = Arrays.copyOf(read.getBaseQualities(), read.getBaseQualities().length);
        for ( int i = startOffset; i < startOffset + length; i++ )
            withLowQuals[i] = (byte)(read.getBaseQualities()[i] + (i % 2 == 0 ? -1 : 0));
        readWithLowQuals.setBaseQualities(withLowQuals);
        return readWithLowQuals;
    }

    private SAMRecord makeOverlappingRead(final String leftFlank, final int leftQual, final String overlapBases,
                                              final byte[] overlapQuals, final String rightFlank, final int rightQual,
                                              final int alignmentStart) {
        final String bases = leftFlank + overlapBases + rightFlank;
        final int readLength = bases.length();
        final SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, alignmentStart, readLength);
        final byte[] leftQuals = Utils.dupBytes((byte) leftQual, leftFlank.length());
        final byte[] rightQuals = Utils.dupBytes((byte) rightQual, rightFlank.length());
        final byte[] quals = Utils.concat(leftQuals, overlapQuals, rightQuals);
        read.setCigarString(readLength + "M");
        read.setReadBases(bases.getBytes());
        read.setBaseQualities(quals);
        ReadUtils.setInsertionBaseQualities(read, quals);
        ReadUtils.setDeletionBaseQualities(read, quals);
        ReadUtils.setReadGroup(read, rgForMerged);
        read.setMappingQuality(60);
        return read;
    }

    @DataProvider(name = "MergeFragmentsOffContig")
    public Object[][] makeMergeFragmentsOffContig() throws Exception {
        List<Object[]> tests = new ArrayList<>();

        for ( final int pre1 : Arrays.asList(0, 50)) {
            for ( final int post1 : Arrays.asList(0, 50)) {
                for ( final int pre2 : Arrays.asList(0, 50)) {
                    for ( final int post2 : Arrays.asList(0, 50)) {
                        tests.add(new Object[]{pre1, post1, pre2, post2});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    private static final byte highQuality = 30;
    private static final byte overlappingQuality = 20;

    @DataProvider(name = "AdjustFragmentsTest")
    public Object[][] createAdjustFragmentsTest() throws Exception {
        List<Object[]> tests = new ArrayList<>();

        final String leftFlank = "CCC";
        final String rightFlank = "AAA";
        final String allOverlappingBases = "ACGTACGTGGAACCTTAG";
        for ( int overlapSize = 1; overlapSize < allOverlappingBases.length(); overlapSize++ ) {
            final String overlappingBases = allOverlappingBases.substring(0, overlapSize);
            final byte[] overlappingBaseQuals = new byte[overlapSize];
            for ( int i = 0; i < overlapSize; i++ ) overlappingBaseQuals[i] = highQuality;
            final SAMRecord read1  = makeOverlappingRead(leftFlank, highQuality, overlappingBases, overlappingBaseQuals, "", highQuality, 1);
            final SAMRecord read2  = makeOverlappingRead("", highQuality, overlappingBases, overlappingBaseQuals, rightFlank, highQuality, leftFlank.length() + 1);
            tests.add(new Object[]{read1, read2, overlapSize});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "AdjustFragmentsTest")
    public void testAdjustingTwoReads(final SAMRecord read1, final SAMRecord read2, final int overlapSize) {
        FragmentUtils.adjustQualsOfOverlappingPairedFragments(Pair.of(read1, read2));
        checkStuff(read1, read2, overlapSize);
    }

    @Test(dataProvider = "AdjustFragmentsTest")
    public void testAdjustingTwoReadsReversed(final SAMRecord read1, final SAMRecord read2, final int overlapSize) {
        FragmentUtils.adjustQualsOfOverlappingPairedFragments(Pair.of(read2, read1));
        checkStuff(read1, read2, overlapSize);
    }

    public void checkStuff(SAMRecord read1, SAMRecord read2, int overlapSize) {
        for ( int i = 0; i < read1.getReadLength() - overlapSize; i++ )
            Assert.assertEquals(read1.getBaseQualities()[i], highQuality);
        for ( int i = read1.getReadLength() - overlapSize; i < read1.getReadLength(); i++ )
            Assert.assertEquals(read1.getBaseQualities()[i], overlappingQuality);

        for ( int i = 0; i < overlapSize; i++ )
            Assert.assertEquals(read2.getBaseQualities()[i], overlappingQuality);
        for ( int i = overlapSize; i < read2.getReadLength(); i++ )
            Assert.assertEquals(read2.getBaseQualities()[i], highQuality);
    }
}
