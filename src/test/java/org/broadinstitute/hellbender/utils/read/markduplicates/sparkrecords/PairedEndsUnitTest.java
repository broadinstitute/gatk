package org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.LibraryIdGenerator;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.markduplicates.util.ReadEnds;

import java.util.Arrays;
import java.util.Collections;

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

        SAMFileHeader header = hg19Header.clone();
        header.setReadGroups(Arrays.asList(new SAMReadGroupRecord("1")));
        primaryRead.setReadGroup("1");
        secondaryRead.setReadGroup("1");


        // Creating a PairedEnds object and asserting the orientation is correct
        Pair pair = PairedEnds.newPair(primaryRead, secondaryRead, header, 0, MarkDuplicatesScoringStrategy.SUM_OF_BASE_QUALITIES, Collections.singletonMap(MarkDuplicatesSparkUtils.getLibraryForRead(primaryRead, header, LibraryIdGenerator.UNKNOWN_LIBRARY), (byte) 0));
        Assert.assertEquals(pair.getOrientationForPCRDuplicates(), PCROrientation);
        Assert.assertEquals(pair.getOrientationForOpticalDuplicates(), opticalOrientation);
    }
}