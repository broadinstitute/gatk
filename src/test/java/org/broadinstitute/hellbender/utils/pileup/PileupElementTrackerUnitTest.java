package org.broadinstitute.hellbender.utils.pileup;

import htsjdk.samtools.SAMUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class PileupElementTrackerUnitTest extends BaseTest {

    @DataProvider(name = "FixPairOverlappingQualitiesTest")
    public static Object[][] overlappingElementsToFix() throws Exception {
        final byte highQual = 60;
        final byte lowQual = 10;
        final byte qualitySum = 70;
        final byte reducedQuality = 48; // 80% of 60
        final byte zeroQuality = 0;
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead("6M");
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead("6M");
        // set the paired and mate state
        read1.setIsPaired(true);
        read1.setMatePosition(read2);
        // set bases and qualities
        read1.setBases(new byte[] {'A', 'A', 'A', 'A', 'T', 'T'});
        read1.setBaseQualities(new byte[] {highQual, highQual, highQual, lowQual,
                highQual, SAMUtils.MAX_PHRED_SCORE});
        // set the paired and mate state
        read2.setIsPaired(true);
        read2.setMatePosition(read1);
        // set bases and qualities
        read2.setBases(new byte[] {'A', 'T', 'T', 'T', 'T', 'T'});
        read2.setBaseQualities(new byte[] {lowQual, highQual, lowQual, highQual,
                highQual, SAMUtils.MAX_PHRED_SCORE});
        return new Object[][] {
                // Same base, first element with higher quality
                {PileupElement.createPileupForReadAndOffset(read1, 0),
                        PileupElement.createPileupForReadAndOffset(read2, 0),
                        qualitySum, zeroQuality},
                // Different base, both with same quality
                {PileupElement.createPileupForReadAndOffset(read1, 1),
                        PileupElement.createPileupForReadAndOffset(read2, 1),
                        reducedQuality, zeroQuality},
                // Different base, first with higher quality
                {PileupElement.createPileupForReadAndOffset(read1, 2),
                        PileupElement.createPileupForReadAndOffset(read2, 2),
                        reducedQuality, zeroQuality},
                // Different base, second with higher quality
                {PileupElement.createPileupForReadAndOffset(read1, 3),
                        PileupElement.createPileupForReadAndOffset(read2, 3),
                        zeroQuality, reducedQuality},
                // Same base, high quality for both (simple cap)
                {PileupElement.createPileupForReadAndOffset(read1, 4),
                        PileupElement.createPileupForReadAndOffset(read2, 4),
                        QualityUtils.MAX_SAM_QUAL_SCORE, zeroQuality},
                // Same base, maximum qualities for both (simple cap)
                {PileupElement.createPileupForReadAndOffset(read1, 5),
                        PileupElement.createPileupForReadAndOffset(read2, 5),
                        QualityUtils.MAX_SAM_QUAL_SCORE, zeroQuality},
        };
    }

    @Test(dataProvider = "FixPairOverlappingQualitiesTest")
    public void testFixPairOverlappingQualities(final PileupElement first,
            final PileupElement second, final byte expectedQualFirst,
            final byte expectedQualSecond) throws Exception {
        PileupElementTracker.fixPairOverlappingQualities(first, second);
        Assert.assertEquals(first.getQual(), expectedQualFirst);
        Assert.assertEquals(second.getQual(), expectedQualSecond);
    }

    @Test
    public void testFixOverlappingQualitiesDeletedElements() throws Exception {
        // Create two reads with deletion
        final GATKRead read = ArtificialReadUtils.createArtificialRead("1M1D1M");
        read.setBases(new byte[] {'A', 'A'});
        read.setBaseQualities(new byte[] {60, 60});
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead("1M1D1M");
        read.setBases(new byte[] {'A', 'A'});
        read2.setBaseQualities(new byte[] {40, 40});
        // Creates PileupElement with a copy of the reads to keep the original read unmodified, because the read from
        // the PileupElement is a reference to the Object passed to the constructor. This is required
        // because fixing the overlapping qualities change the read stored into PileupElement and we have to check that
        // it does not change anything from the read.
        final GATKRead copy = read.deepCopy();
        final GATKRead copy2 = read2.deepCopy();
        final PileupElement deleted = new PileupElement(copy, 1, copy.getCigarElement(1), 1, 0);
        final PileupElement normal = new PileupElement(copy2, 0, copy.getCigarElement(0), 0, 0);
        PileupElementTracker.fixPairOverlappingQualities(deleted, normal);
        // Copy reads are a reference to the read into the PileupElement, and thus should not be modified.
        // Thus, we use them to test if it is not modified. This could be checked uncommenting the next two lines:
        // Assert.assertSame(deleted.getRead(), copy);
        // Assert.assertSame(normal.getRead(), copy2);
        Assert.assertEquals(copy, read);
        Assert.assertEquals(copy2, read2);
        PileupElementTracker.fixPairOverlappingQualities(normal, deleted);
        Assert.assertEquals(copy, read);
        Assert.assertEquals(copy2, read2);
    }

    @Test
    public void testFixPairOverlappingQualitiesCap() {
        final PileupElement element1 = PileupElement
                .createPileupForReadAndOffset(ArtificialReadUtils.createArtificialRead("1M"), 0);
        element1.getRead().setName("read1");
        element1.getRead().setBases(new byte[] {'A'});
        final PileupElement element2 = PileupElement
                .createPileupForReadAndOffset(ArtificialReadUtils.createArtificialRead("1M"), 0);
        element2.getRead().setName("read2");
        element2.getRead().setBases(new byte[] {'A'});
        final PileupElement element3 = PileupElement
                .createPileupForReadAndOffset(ArtificialReadUtils.createArtificialRead("1M"), 0);
        element3.getRead().setName("read3");
        element2.getRead().setBases(new byte[] {'A'});
        element3.getRead().setBaseQualities(new byte[] {Byte.MAX_VALUE});
        // Check all possible combinations that goes beyond the maximum value
        for (byte i = 0; i < Byte.MAX_VALUE; i++) {
            final byte[] iArray = new byte[] {i};
            final byte j = (byte) (Byte.MAX_VALUE - i);
            element1.getRead().setBaseQualities(iArray);
            element2.getRead().setBaseQualities(new byte[] {j});
            logger.debug("Test: fixing ({}) and ({})", element1, element2);
            PileupElementTracker.fixPairOverlappingQualities(element1, element2);
            Assert.assertEquals(element1.getQual(), QualityUtils.MAX_SAM_QUAL_SCORE);
            Assert.assertEquals(element2.getQual(), 0);
            element1.getRead().setBaseQualities(iArray);
            logger.debug("Test: fixing ({}) and ({})", element3, element1);
            PileupElementTracker.fixPairOverlappingQualities(element3, element1);
            Assert.assertEquals(element3.getQual(), QualityUtils.MAX_SAM_QUAL_SCORE);
            Assert.assertEquals(element1.getQual(), 0);
        }
        // Finally, check what happens if both are the maximum value
        element1.getRead().setBaseQualities(new byte[] {Byte.MAX_VALUE});
        logger.debug("Test: fixing ({}) and ({})", element3, element1);
        PileupElementTracker.fixPairOverlappingQualities(element3, element1);
        Assert.assertEquals(element3.getQual(), QualityUtils.MAX_SAM_QUAL_SCORE);
        Assert.assertEquals(element1.getQual(), 0);
    }
}