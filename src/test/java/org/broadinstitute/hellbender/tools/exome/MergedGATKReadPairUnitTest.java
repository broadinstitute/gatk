package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.*;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.RandomDNA;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Unit tests for {@link MergedGATKReadPair}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class MergedGATKReadPairUnitTest {

    // length of each test read.
    private final int READ_LENGTH = 100;

    // average length of the fragment
    // must be large enough so that some reads pairs wont overlap
    // but small enough so that some will.
    private final int FRAGMENT_LENGTH_AVG = 200;

    // sd of the fragment length
    private final double FRAGMENT_LENGTH_SD = 50;

    // name of the reference chromosome.
    private final String REFERENCE_CHR = "1";

    // probability of a clip at either end of a read.
    private final double CLIP_PROB = 0.1;

    // determine the rate of a clip to finish (determines the length of the clips).
    private final double CLIP_TO_NON_CLIP_PROB = 0.05;

    // transition prob from clip/ins/del to matching block:
    private final double INDEL_TO_MATCH_PROB = 0.5;

    // transition prob from matching to ins/del
    private final double MATCH_TO_INDEL_PROB = 0.01;

    // base call error prob.
    private final double ERROR_PROB = 0.01;

    private final int TEST_REPEATS = 100;

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = IllegalArgumentException.class)
    public void testUnmappedRight(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left, final GATKRead right) {
        final GATKRead unmappedRight = right.copy();
        unmappedRight.setIsUnmapped();
        MergedGATKReadPair.mergeReadPair(left, unmappedRight);
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = IllegalArgumentException.class)
    public void testUnmappedLeft(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                 final GATKRead right) {
        final GATKRead unmappedLeft = left.copy();
        unmappedLeft.setIsUnmapped();
        MergedGATKReadPair.mergeReadPair(unmappedLeft, right);
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testSwappedLeftAndRight(final byte[] fragment, final GATKRead left, final GATKRead right) {
        testMergeOutcome(fragment, right, left);
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testNonSwappedLeftAndRight(final byte[] fragment, final GATKRead left, final GATKRead right) {
        testMergeOutcome(fragment, left, right);
    }

    public void testMergeOutcome(final byte[] fragment, final GATKRead left, final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertSame(left.getName(), merged.getName());
        Assert.assertSame(left.getReadGroup(), merged.getReadGroup());
        Assert.assertEquals(left.getContig(), merged.getContig());
        Assert.assertEquals(Math.min(left.getStart(), right.getStart()), merged.getStart());
        Assert.assertEquals(merged.getUnclippedStart(), merged.getStart());
        Assert.assertTrue(merged.getCigar().getCigarElements().size() > 0);
        Assert.assertFalse(merged.getCigar().getCigarElements().stream().anyMatch(ce -> ce.getOperator() != CigarOperator.D && ce.getOperator() != CigarOperator.M));
        for (int i = 1; i < merged.getCigar().getCigarElements().size(); i++) {
            Assert.assertNotEquals(merged.getCigar().getCigarElements().get(i).getOperator(),
                    merged.getCigar().getCigarElements().get(i - 1).getOperator());
        }

        // read-pair to reference base and quality score projection.
        // the next two for loop bellow will project the read base-calls into this arrays in reference coordinate
        // space. The last loop check that the merged read in consistent with the projection.
        final byte[] bases = new byte[fragment.length * 30]; // 30 large number  to make sure that the reference fits in
        final byte[] quals = new byte[fragment.length * 30];
        // Filling the read-pair call to reference projection with -1 indicating that neither read aligns
        // to that position; by default and before we start to mapping the reads onto the reference that is the case
        // at every position.
        Arrays.fill(bases,(byte) -1);


        int referenceNextPosition = left.getStart();
        int readNextPosition = 0;
        // copy the left read calls as they are.
        for (final CigarElement element : left.getCigar().getCigarElements()) {
            switch (element.getOperator()) {
                case M:
                case X:
                case EQ:
                    System.arraycopy(left.getBases(), readNextPosition, bases, referenceNextPosition, element.getLength());
                    System.arraycopy(left.getBaseQualities(), readNextPosition, quals, referenceNextPosition, element.getLength());
                    readNextPosition += element.getLength();
                    referenceNextPosition += element.getLength();
                    break;
                case I:
                case S:
                case H:
                    readNextPosition += element.getLength();
                    break;
                case D:
                    referenceNextPosition += element.getLength();
                    break;
            }
        }

        // overlay right read calls considering mismatches.
        referenceNextPosition = right.getStart();
        readNextPosition = 0;
        for (final CigarElement element : right.getCigar().getCigarElements()) {
            switch (element.getOperator()) {
                case M:
                case X:
                case EQ:
                    final byte[] rightBases = Arrays.copyOfRange(right.getBases(), readNextPosition, readNextPosition + element.getLength());
                    final byte[] rightQuals = Arrays.copyOfRange(right.getBaseQualities(), readNextPosition, readNextPosition + element.getLength());
                    for (int i = 0; i < element.getLength(); i++) {
                        if (rightQuals[i] > quals[i + referenceNextPosition]) {
                            quals[i + referenceNextPosition] = rightQuals[i];
                            bases[i + referenceNextPosition] = rightBases[i];
                        }
                    }
                    readNextPosition += element.getLength();
                    referenceNextPosition += element.getLength();
                    break;
                case I:
                case S:
                case H:
                    readNextPosition += element.getLength();
                    break;
                case D:
                    referenceNextPosition += element.getLength();
                    break;
            }
        }

        /// Then we compare the resulting merged read with the reference bases and quals read-pair projection.
        referenceNextPosition = merged.getStart();
        readNextPosition = 0;
        for (final byte prefix : Arrays.copyOfRange(bases, 0, referenceNextPosition)) {
            Assert.assertEquals(prefix, (byte) -1);
        }
        for (final CigarElement element : merged.getCigar().getCigarElements()) {
            switch (element.getOperator()) {
                case M:
                    Assert.assertEquals(
                            Arrays.copyOfRange(merged.getBases(), readNextPosition, readNextPosition + element.getLength()),
                            Arrays.copyOfRange(bases, referenceNextPosition, referenceNextPosition + element.getLength()));
                    readNextPosition += element.getLength();
                    referenceNextPosition += element.getLength();
                    break;
                case D:
                    for (int i = 0; i < element.getLength(); i++) {
                        Assert.assertEquals(bases[referenceNextPosition + i], (byte) -1);
                    }
                    referenceNextPosition += element.getLength();
            }
        }

        for (int i = referenceNextPosition; i < bases.length; i++) {
            Assert.assertEquals(bases[i], (byte) - 1);
        }
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testName(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left, final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertEquals(merged.getName(), left.getName());
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetName(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left, final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setName("BLAH");
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testLength(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left, final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertEquals(merged.getLength(), merged.getBases().length);
    }


    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetPosition(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left, final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setPosition(merged.getContig(), merged.getStart() + 1);
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetPositionLocatable(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left, final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setPosition(new SimpleInterval(merged.getContig(), merged.getStart(), merged.getStart()));
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testGetUnclippedStart(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                      final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertEquals(merged.getUnclippedStart(), left.getStart());
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testGetUnclippedEnd(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                      final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertEquals(merged.getUnclippedEnd(), right.getEnd());
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testGetMateContig(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                    final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertNull(merged.getMateContig());
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testGetMateStart(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                  final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertEquals(merged.getMateStart(), 0);
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetMatePosition(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                 final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setMatePosition(merged.getContig(), 1);
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetMatePositionLocatable(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                    final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setMatePosition(new SimpleInterval(merged.getContig(), 1, 1));
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testGetFragmentLength(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                   final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertEquals(merged.getFragmentLength(), 0);
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetFragmentLength(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                             final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setFragmentLength(100);
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testGetMappingQuality(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                   final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertEquals(merged.getMappingQuality(), Math.max(left.getMappingQuality(), right.getMappingQuality()));
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetMappingQuality(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                      final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setMappingQuality(100);
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testGetBases(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                   final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertEquals(merged.getBases().length, merged.getLength());
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetBases(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                      final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setBases(new byte[merged.getBases().length]);
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testGetBaseQuality(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                     final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertEquals(merged.getBaseQualities().length, merged.getLength());
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetBaseQualities(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                             final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setBaseQualities(new byte[merged.getBases().length]);
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testGetCigar(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                     final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertNotNull(merged.getCigar());
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetCigar(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                     final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setCigar(merged.getCigar());
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testGetReadGroup(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                   final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertEquals(merged.getReadGroup(), null);
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetReadGroup(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                             final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setReadGroup("OTHER_GROUP");
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testIsPaired(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                     final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertFalse(merged.isPaired());
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetIsPaired(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                 final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setIsPaired(true);
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testIsProperlyPaired(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                   final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertFalse(merged.isProperlyPaired());
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetIsProperlyPaired(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setIsProperlyPaired(true);
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testIsUnmapped(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                   final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertFalse(merged.isUnmapped());
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetIsUnmapped(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                        final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setIsUnmapped();
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testMateIsUnmapped(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                    final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertFalse(merged.mateIsUnmapped());
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetMateIsUnmapped(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                  final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setMateIsUnmapped();
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testIsReverseStrand(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                        final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertFalse(merged.isReverseStrand());
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetIsReverseStrand(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                      final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setIsReverseStrand(true);
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testMateIsReverseStrand(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                   final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertFalse(merged.mateIsReverseStrand());
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetMateIsReverseStrand(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                       final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setMateIsReverseStrand(true);
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testIsFirstInPair(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                   final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertFalse(merged.isFirstOfPair());
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetIsFirstOfPair(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                       final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setIsFirstOfPair();
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testIsSecondInPair(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                             final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertFalse(merged.isSecondOfPair());
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetIsSecondOfPair(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                     final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setIsSecondOfPair();
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testIsSecondaryAlignment(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                             final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertFalse(merged.isSecondaryAlignment());
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetIsSecondaryAlignment(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                     final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setIsSecondaryAlignment(true);
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testIsSupplementaryAlignment(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                 final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertFalse(merged.isSupplementaryAlignment());
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetIsSupplementaryAlignment(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                     final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setIsSupplementaryAlignment(true);
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testFailsVendor(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                 final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertFalse(merged.failsVendorQualityCheck());
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetFailsVendor(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                     final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setFailsVendorQualityCheck(true);
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testIsDuplicated(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                    final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertFalse(merged.isDuplicate());
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetIsDuplicated(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                     final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setIsDuplicate(true);
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testHasAttribute(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                     final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.hasAttribute("PROB");
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testGetAttributeAsInteger(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                     final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.getAttributeAsInteger("PROB");
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testGetAttributeAsString(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                          final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.getAttributeAsString("PROB");
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testGetAttributeAsByteArray(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                         final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.getAttributeAsByteArray("PROB");
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetAttributeInteger(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                         final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setAttribute("PROB", 1);
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetAttributeString(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                        final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setAttribute("PROB", "VALUE");
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testSetAttributeByteArray(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                       final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.setAttribute("PROB", new byte[merged.getLength()]);
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testClearAttribute(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                          final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.clearAttribute("PROB");
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testClearAttributes(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                          final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.clearAttributes();
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testCopy(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                             final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        final GATKRead copy = merged.copy();
        Assert.assertEquals(copy.getStart(), merged.getStart());
        Assert.assertSame(copy.getBases(), merged.getBases());
        Assert.assertSame(copy.getBaseQualities(), merged.getBaseQualities());
        Assert.assertEquals(copy.getLength(), merged.getLength());
        Assert.assertSame(copy.getName(), merged.getName());
        Assert.assertSame(copy.getCigar(), merged.getCigar());
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testDeepCopy(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                        final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        final GATKRead copy = merged.deepCopy();
        Assert.assertEquals(copy.getStart(), merged.getStart());
        Assert.assertEquals(copy.getBases(), merged.getBases());
        Assert.assertEquals(copy.getBaseQualities(), merged.getBaseQualities());
        Assert.assertEquals(copy.getLength(), merged.getLength());
        Assert.assertEquals(copy.getName(), merged.getName());
        Assert.assertEquals(copy.getCigar(), merged.getCigar());
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testConvertToSAMRecord(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                          final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.convertToSAMRecord(new SAMFileHeader());
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testConvertToGoogle(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                       final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.convertToGoogleGenomicsRead();
    }

    @Test(dataProvider = "randomFragmentAndReadsData", expectedExceptions = UnsupportedOperationException.class)
    public void testGetSAMString(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                    final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        merged.getSAMString();
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testContig(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                                 final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertEquals(merged.getContig(), REFERENCE_CHR);
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testStart(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                           final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertEquals(merged.getStart(), left.getStart());
    }

    @Test(dataProvider = "randomFragmentAndReadsData")
    public void testEnd(@SuppressWarnings("unused") final byte[] fragment, final GATKRead left,
                          final GATKRead right) {
        final GATKRead merged = MergedGATKReadPair.mergeReadPair(left, right);
        Assert.assertEquals(merged.getEnd(), right.getEnd());
    }

    @DataProvider(name = "randomFragmentAndReadsData")
    private Object[][] randomFragmentAndReadsData() {
        final SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord(REFERENCE_CHR, 3000));
        final Random rdn = new Random(1313);
        final List<Object[]> result = new ArrayList<>(TEST_REPEATS);
        for (int i = 0; i < TEST_REPEATS; i++) {
            final int fragmentSize = (int) Math.max(READ_LENGTH, rdn.nextGaussian() * FRAGMENT_LENGTH_SD + FRAGMENT_LENGTH_AVG);
            final byte[] fragment = new RandomDNA(rdn).nextBases(fragmentSize);
            final Pair<GATKRead, GATKRead> readPair = randomReadPair(fragment, rdn, header);
            result.add(new Object[] { fragment, readPair.getLeft(), readPair.getRight() });
        }
        return result.toArray(new Object[result.size()][]);
    }

    private Pair<GATKRead, GATKRead> randomReadPair(final byte[] fragment, final Random rdn, final SAMFileHeader header) {
        final Cigar fragmentCigar = randomCigar(fragment.length, rdn);
        final Cigar leftCigar = leftCigar(fragmentCigar);
        final int leftStart = leftStart(leftCigar);
        final Cigar rightCigar = rightCigar(fragmentCigar);
        final int rightStart = rightStart(leftStart, fragmentCigar, rightCigar);
        final SAMRecord leftRead = new SAMRecord(header);
        final String name = "READ_PAIR_" + rdn.nextInt(TEST_REPEATS * 100);
        leftRead.setReadName(name);
        leftRead.setReadBases(introduceErrors(Arrays.copyOfRange(fragment, 0, READ_LENGTH), rdn));
        leftRead.setBaseQualities(randomQuals(READ_LENGTH, rdn));
        leftRead.setCigar(leftCigar);
        leftRead.setReferenceIndex(0);
        leftRead.setReferenceName(REFERENCE_CHR);
        leftRead.setMappingQuality(rdn.nextInt(60) + 1);
        leftRead.setAlignmentStart(leftStart);
        final GATKRead leftGATKRead = new SAMRecordToGATKReadAdapter(leftRead);
        final SAMRecord rightRead = new SAMRecord(header);
        rightRead.setReadBases(introduceErrors(Arrays.copyOfRange(fragment, fragment.length - READ_LENGTH, fragment.length), rdn));
        rightRead.setReadName(name);
        rightRead.setMappingQuality(rdn.nextInt(60) + 1);
        rightRead.setCigar(rightCigar);
        rightRead.setBaseQualities(randomQuals(READ_LENGTH, rdn));
        rightRead.setReferenceIndex(0);
        rightRead.setReferenceName(REFERENCE_CHR);
        rightRead.setAlignmentStart(rightStart);
        final GATKRead rightGATKRead = new SAMRecordToGATKReadAdapter(rightRead);
        return new ImmutablePair<>(leftGATKRead, rightGATKRead);
    }

    private byte[] randomQuals(final int readLength, final Random rdn) {
        final byte[] quals = new byte[readLength];
        for (int i = 0; i < readLength; i++) {
            quals[i] = (byte) (10 + rdn.nextInt(20));
        }
        return quals;
    }

    private byte[] introduceErrors(final byte[] bases, final Random rdn) {
        for (int i = 0; i < bases.length; i++) {
            if (rdn.nextDouble() < ERROR_PROB) {
                bases[i] = new RandomDNA(rdn).nextBases(1)[0];
            }
        }
        return bases;
    }

    private int leftStart(final Cigar leftCigar) {
        int referenceStart = 1;
        for (final CigarElement element : leftCigar.getCigarElements()) {
            if (element.getOperator().consumesReadBases() && element.getOperator().consumesReferenceBases()) {
                return referenceStart;
            } else if (element.getOperator().consumesReferenceBases()) {
                referenceStart += element.getLength();
            }
        }
        return 1;
    }

    private int rightStart(final int leftStart, final Cigar fragmentCigar, final Cigar rightCigar) {
        final int fragmentReferenceLength = referenceLength(fragmentCigar);
        final int rightReferenceLength = referenceLength(rightCigar);
        return leftStart + fragmentReferenceLength - rightReferenceLength;
    }

    private int referenceLength(final Cigar cigar) {
        int result = 0;
        for (final CigarElement element : cigar.getCigarElements()) {
            if (element.getOperator().consumesReferenceBases()) {
                result += element.getLength();
            }
        }
        return result;
    }

    private Cigar leftCigar(final Cigar fragmentCigar) {
        return new Cigar(subsetReadCigarElements(fragmentCigar.getCigarElements()));
    }

    private Cigar rightCigar(final Cigar fragmentCigar) {
        final List<CigarElement> fragmentCigarElements = new ArrayList<>(fragmentCigar.getCigarElements());
        Collections.reverse(fragmentCigarElements);
        final List<CigarElement> result = subsetReadCigarElements(fragmentCigarElements);
        Collections.reverse(result);
        return new Cigar(result);
    }

    private List<CigarElement> subsetReadCigarElements(final List<CigarElement> input) {
        int remainingReadLength = READ_LENGTH;
        final List<CigarElement> result = new ArrayList<>(input.size());
        for (final CigarElement next : input) {
            if (!next.getOperator().consumesReadBases()) {
                result.add(next);
            } else if (next.getLength() >= remainingReadLength) {
                result.add(new CigarElement(remainingReadLength, next.getOperator()));
                remainingReadLength = 0;
            } else {
                result.add(new CigarElement(next.getLength(), next.getOperator()));
                remainingReadLength -= next.getLength();
            }
            if (remainingReadLength == 0) {
                break;
            }
        }
        if (result.get(result.size() - 1).getOperator() == CigarOperator.I) {
            result.set(result.size() - 1, new CigarElement(result.get(result.size() - 1).getLength(), CigarOperator.S));
        }
        return result;
    }

    private Cigar randomCigar(final int fragmentLength, final Random rdn) {
        final CigarElement leftClip  = randomClip(rdn);
        final CigarElement rightClip = randomClip(rdn);
        int nonClipLength = fragmentLength - (leftClip == null ? 0 : leftClip.getLength())
                - (rightClip == null ? 0 : rightClip.getLength());
        if (nonClipLength < 1) { // default to all M cigar, should happen rarely.
            return new Cigar(Collections.singletonList(new CigarElement(fragmentLength, CigarOperator.M)));
        }
        final List<CigarElement> elements = new ArrayList<>(10);
        if (leftClip != null) {
            elements.add(leftClip);
        }
        CigarOperator operator = CigarOperator.M;
        int length = 0;
        while (nonClipLength > 0) {
            length++;
            if (operator.consumesReadBases()) {
                nonClipLength--;
            }
            switch (operator) {
                case M:
                    if (rdn.nextDouble() <= MATCH_TO_INDEL_PROB) {
                        elements.add(new CigarElement(length, CigarOperator.M));
                        operator = rdn.nextDouble() < .5 ? CigarOperator.D : CigarOperator.I;
                        length = 0;
                    }
                    break;
                case I:
                    if (rdn.nextDouble() <= INDEL_TO_MATCH_PROB) {
                        elements.add(new CigarElement(length, CigarOperator.I));
                        operator = CigarOperator.M;
                        length = 0;
                    }
                    break;
                case D:
                    if (rdn.nextDouble() <= INDEL_TO_MATCH_PROB) {
                        elements.add(new CigarElement(length, CigarOperator.D));
                        operator = CigarOperator.M;
                        length = 0;
                    }
            }
        }
        elements.add(new CigarElement(length, operator));
        if (rightClip != null) {
            elements.add(rightClip);
        }
        return new Cigar(elements);
    }

    private CigarElement randomClip(final Random rdn) {
           if (rdn.nextDouble() > CLIP_PROB) {
               return null;
           } else {
               int length = 1;
               while (rdn.nextDouble() > CLIP_TO_NON_CLIP_PROB) {
                   length++;
               }
               return new CigarElement(length, CigarOperator.SOFT_CLIP);
           }
    }
}
