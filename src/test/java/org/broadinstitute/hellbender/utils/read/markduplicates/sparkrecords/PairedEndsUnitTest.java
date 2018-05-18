package org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadEnds;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class PairedEndsUnitTest extends GATKBaseTest {

    @DataProvider
    public Object[][] orientationTruthTable() {
        return new Object[][]{
                {false, false, false, ReadEnds.FF, ReadEnds.FF},
                {false, false, true, ReadEnds.FR, ReadEnds.FR},
                {false, true, false, ReadEnds.RF, ReadEnds.RF},
                {false, true, true, ReadEnds.RR, ReadEnds.RR},
                {true, false, false, ReadEnds.FF, ReadEnds.FF},
                {true, false, true, ReadEnds.RF, ReadEnds.FR},
                {true, true, false, ReadEnds.FR, ReadEnds.RF},
                {true, true, true, ReadEnds.RR, ReadEnds.RR},};
    }

    @Test (dataProvider = "orientationTruthTable")
    public void testGetOrientationForPCRDuplicates(boolean flipStarts, boolean firstReadReverse, boolean secondReadReverse,
                                                   byte PCROrientation, byte opticalOrientation) {
        GATKRead primaryRead;
        GATKRead secondaryRead;

        if (flipStarts) {
            secondaryRead = ArtificialReadUtils.createSamBackedRead("100M",100000,100);
            primaryRead = ArtificialReadUtils.createSamBackedRead("100M",101000,100);
        } else {
            primaryRead = ArtificialReadUtils.createSamBackedRead("100M",100000,100);
            secondaryRead = ArtificialReadUtils.createSamBackedRead("100M",101000,100);
        }
        primaryRead.setIsFirstOfPair();
        primaryRead.setIsReverseStrand(firstReadReverse);
        secondaryRead.setIsSecondOfPair();
        secondaryRead.setIsReverseStrand(secondReadReverse);


        // Creating a PairedEnds object and asserting the orientation is correct
        Pair pair = PairedEnds.newPair(primaryRead, secondaryRead, hg19Header, 0, MarkDuplicatesScoringStrategy.SUM_OF_BASE_QUALITIES);
        Assert.assertEquals(pair.getOrientationForPCRDuplicates(), PCROrientation);
        Assert.assertEquals(pair.getOrientationForOpticalDuplicates(), opticalOrientation);
    }
}