package org.broadinstitute.hellbender.utils.fragments;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class FragmentUtilsUnitTest extends GATKBaseTest {

    private static final byte HIGH_QUALITY = 30;
    private static final byte OVERLAPPING_QUALITY = 20;

    private GATKRead makeOverlappingRead(final String leftFlank, final int leftQual, final String overlapBases,
                                              final byte[] overlapQuals, final String rightFlank, final int rightQual,
                                              final int alignmentStart, final int leftSoftclip, final int rightSoftclip) {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        header.addReadGroup(new SAMReadGroupRecord("RG1"));
        final String bases = leftFlank + overlapBases + rightFlank;
        final int readLength = bases.length();
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, alignmentStart + leftSoftclip, readLength);
        final byte[] leftQuals = Utils.dupBytes((byte) leftQual, leftFlank.length());
        final byte[] rightQuals = Utils.dupBytes((byte) rightQual, rightFlank.length());
        final byte[] quals = Utils.concat(leftQuals, overlapQuals, rightQuals);
        read.setCigar((leftSoftclip != 0 ? leftSoftclip + "S" : "") + (readLength - rightSoftclip - leftSoftclip) + "M" + (rightSoftclip != 0 ? rightSoftclip + "S" : "") );
        read.setBases(bases.getBytes());
        read.setBaseQualities(quals);
        read.setReadGroup("RG1");
        read.setMappingQuality(60);
        return read;
    }

    @DataProvider(name = "AdjustFragmentsTest")
    public Object[][] createAdjustFragmentsTest() throws Exception {
        List<Object[]> tests = new ArrayList<>();

        final String leftFlank = "CCC";
        final String rightFlank = "AAA";
        final String allOverlappingBases = "ACGTACGTGGAACCTTAG";
        for ( int overlapSize = 1; overlapSize < allOverlappingBases.length(); overlapSize++ ) {
            final String overlappingBases = allOverlappingBases.substring(0, overlapSize);
            final byte[] overlappingBaseQuals = new byte[overlapSize];
            for ( int i = 0; i < overlapSize; i++ ) {
                overlappingBaseQuals[i] = HIGH_QUALITY;
            }
            final GATKRead read1  = makeOverlappingRead(leftFlank, HIGH_QUALITY, overlappingBases, overlappingBaseQuals, "", HIGH_QUALITY, 1, 0, 0 );
            final GATKRead read2  = makeOverlappingRead("", HIGH_QUALITY, overlappingBases, overlappingBaseQuals, rightFlank, HIGH_QUALITY, leftFlank.length() + 1, 0, 0);
            tests.add(new Object[]{read1, read2, overlapSize});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "AdjustFragmentsTest")
    public void testAdjustingTwoReads(final GATKRead read1, final GATKRead read2, final int overlapSize) {
        FragmentUtils.adjustQualsOfOverlappingPairedFragments(ImmutablePair.of(read1, read2), true, OptionalInt.empty(), OptionalInt.empty());

        for ( int i = 0; i < read1.getLength() - overlapSize; i++ ) {
            Assert.assertEquals(read1.getBaseQualities()[i], HIGH_QUALITY);
        }
        for ( int i = read1.getLength() - overlapSize; i < read1.getLength(); i++ ) {
            Assert.assertEquals(read1.getBaseQualities()[i], OVERLAPPING_QUALITY);
        }

        for ( int i = 0; i < overlapSize; i++ ) {
            Assert.assertEquals(read2.getBaseQualities()[i], OVERLAPPING_QUALITY);
        }
        for ( int i = overlapSize; i < read2.getLength(); i++ ) {
            Assert.assertEquals(read2.getBaseQualities()[i], HIGH_QUALITY);
        }
    }


    // Generate a bunch of reads with softclips that do not overlap with the other read.
    @DataProvider(name = "AdjustFragmentsTestSoftClipsNotOverlapping")
    public Object[][] createAdjustFragmentsTestSoftClips() throws Exception {
        List<Object[]> tests = new ArrayList<>();

        final String leftFlank = "CCC";
        final String rightFlank = "AAA";
        final String allOverlappingBases = "ACGTACGTGGAACCTTAG";
        for ( int overlapSize = 1; overlapSize < allOverlappingBases.length(); overlapSize++ ) {
            for (int leftSoftclip = 0; leftSoftclip <= leftFlank.length(); leftSoftclip++) {
                for (int rightSoftclip = 0; rightSoftclip <= rightFlank.length(); rightSoftclip++) {
                    final String overlappingBases = allOverlappingBases.substring(0, overlapSize);
                    final byte[] overlappingBaseQuals = new byte[overlapSize];
                    for ( int i = 0; i < overlapSize; i++ ) {
                        overlappingBaseQuals[i] = HIGH_QUALITY;
                    }
                    final GATKRead read1  = makeOverlappingRead(leftFlank, HIGH_QUALITY, overlappingBases, overlappingBaseQuals, "", HIGH_QUALITY, 1, leftSoftclip, 0);
                    final GATKRead read2  = makeOverlappingRead("", HIGH_QUALITY, overlappingBases, overlappingBaseQuals, rightFlank, HIGH_QUALITY, leftFlank.length() + 1, 0, rightSoftclip);
                    tests.add(new Object[]{read1, read2, overlapSize});
                }
            }
        }
        return tests.toArray(new Object[][]{});
    }

    // Assert that despite the softclips that the reads are being properly
    @Test(dataProvider = "AdjustFragmentsTestSoftClipsNotOverlapping")
    public void testAdjustingTwoReadsWithSoftClipping(final GATKRead read1, final GATKRead read2, final int overlapSize) {
        FragmentUtils.adjustQualsOfOverlappingPairedFragments(ImmutablePair.of(read1, read2), true, OptionalInt.empty(), OptionalInt.empty());

        for ( int i = 0; i < read1.getLength() - overlapSize; i++ ) {
            Assert.assertEquals(read1.getBaseQualities()[i], HIGH_QUALITY);
        }
        for ( int i = read1.getLength() - overlapSize; i < read1.getLength(); i++ ) {
            Assert.assertEquals(read1.getBaseQualities()[i], OVERLAPPING_QUALITY);
        }

        for ( int i = 0; i < overlapSize; i++ ) {
            Assert.assertEquals(read2.getBaseQualities()[i], OVERLAPPING_QUALITY);
        }
        for ( int i = overlapSize; i < read2.getLength(); i++ ) {
            Assert.assertEquals(read2.getBaseQualities()[i], HIGH_QUALITY);
        }
    }


    // Generate a bunch of reads with softclips that are overlapping with the other read (and allow for reads with no overlap at all except for softclips)
    @DataProvider(name = "AdjustFragmentsTestSoftClipsInOverlapRegion")
    public Object[][] createAdjustFragmentsTestSoftClipsInOverlapRegion() throws Exception {
        List<Object[]> tests = new ArrayList<>();

        final String leftFlank = "CCC";
        final String rightFlank = "AAA";
        final String allOverlappingBases = "ACGTACGTGGAACCTTAG";
        for ( int overlapSize = 1; overlapSize < allOverlappingBases.length(); overlapSize++ ) {
            for (int leftSoftclip = 0; leftSoftclip <= overlapSize; leftSoftclip++) {
                for (int rightSoftclip = 0; rightSoftclip <= overlapSize; rightSoftclip++) {
                    final String overlappingBases = allOverlappingBases.substring(0, overlapSize);
                    final byte[] overlappingBaseQuals = new byte[overlapSize];
                    for ( int i = 0; i < overlapSize; i++ ) {
                        overlappingBaseQuals[i] = HIGH_QUALITY;
                    }
                    // Flipped so that the softclips occur in the overlapping region instead
                    final GATKRead read1  = makeOverlappingRead(leftFlank, HIGH_QUALITY, overlappingBases, overlappingBaseQuals, "", HIGH_QUALITY, 1, 0, rightSoftclip);
                    final GATKRead read2  = makeOverlappingRead("", HIGH_QUALITY, overlappingBases, overlappingBaseQuals, rightFlank, HIGH_QUALITY, leftFlank.length() + 1, leftSoftclip, 0);
                    tests.add(new Object[]{read1, rightSoftclip, read2, leftSoftclip, overlapSize});
                }
            }
        }
        return tests.toArray(new Object[][]{});
    }

    // Assert that despite the softclips that the reads are being properly
    @Test(dataProvider = "AdjustFragmentsTestSoftClipsInOverlapRegion")
    public void testAdjustingTwoReadsWithSoftClippingOverlapingEachother(final GATKRead read1, final int read1Softclips, final GATKRead read2, final int read2Softclips, final int overlapSize) {
        FragmentUtils.adjustQualsOfOverlappingPairedFragments(ImmutablePair.of(read1, read2), true, OptionalInt.empty(), OptionalInt.empty());

        // Untouched bases in R1 that don't overlap
        for ( int i = 0; i < read1.getLength() - overlapSize; i++ ) {
            Assert.assertEquals(read1.getBaseQualities()[i], HIGH_QUALITY);
        }
        // Bases in R1 that overlap with R2 (before R1 softclips and sans R2 softclips)
        for ( int i = read1.getLength() - overlapSize + read2Softclips; i < read1.getLength() - read1Softclips; i++ ) {
            Assert.assertEquals(read1.getBaseQualities()[i], OVERLAPPING_QUALITY);
        }
        // Softclipped bases from R1 that did not get adjusted
        for ( int i = read1.getLength() - read1Softclips ; i < read1.getLength(); i++ ) {
            Assert.assertEquals(read1.getBaseQualities()[i], HIGH_QUALITY);
        }


        // Softclipped bases from R2 that did not get adjusted
        for ( int i = 0; i < read2Softclips; i++ ) {
            Assert.assertEquals(read2.getBaseQualities()[i], HIGH_QUALITY);
        }
        // Bases in R2 that overlap with R1 (after R2 softclips adjusted for R1 softclips)
        for ( int i = read2Softclips; i < overlapSize - read1Softclips; i++ ) {
            Assert.assertEquals(read2.getBaseQualities()[i], OVERLAPPING_QUALITY);
        }
        // Bases in R2 that did not overlap at all with R1
        for ( int i = overlapSize; i < read2.getLength(); i++ ) {
            Assert.assertEquals(read2.getBaseQualities()[i], HIGH_QUALITY);
        }
    }
}
