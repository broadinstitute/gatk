package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class SVIntervalUtilsUnitTest extends GATKBaseTest {

    @Test(groups = "sv")
    public void getPaddedIntervalTest() {
        final SVInterval interval = new SVInterval(0, 100, 200);
        final int padding = 10;
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("seq0", 1000));
        final SVInterval paddedInterval = SVIntervalUtils.getPaddedInterval(interval, padding, dictionary);
        Assert.assertEquals(paddedInterval.getContig(), 0);
        Assert.assertEquals(paddedInterval.getStart(), 90);
        Assert.assertEquals(paddedInterval.getEnd(), 210);

        final SVInterval interval2 = new SVInterval(0, 5, 995);
        final SVInterval paddedInterval2 = SVIntervalUtils.getPaddedInterval(interval2, padding, dictionary);
        Assert.assertEquals(paddedInterval2.getContig(), 0);
        Assert.assertEquals(paddedInterval2.getStart(), 0);
        Assert.assertEquals(paddedInterval2.getEnd(), 999);
    }

    @DataProvider(name = "reciprocalOverlapTest")
    public Object[][] getReciprocalOverlapData() {
        return new Object[][]{
                {10, 20, 0, 9, 0., 1.},
                {10, 20, 20, 29, 0., 1.},
                {10, 20, 5, 15, 5. / 10., 6. / 10.},
                {10, 20, 15, 25, 5. / 10., 6. / 10.},
                {10, 20, 19, 20, 1. / 10., 2. / 10.},
                {10, 20, 19, 21, 1. / 10., 2. / 10.},
                {10, 20, 18, 21, 2. / 10., 3. / 10.},
                {10, 20, 18, 30, 2. / 12., 3. / 12.}
        };
    }

    @Test(groups = "sv",
            dataProvider = "reciprocalOverlapTest")
    public void reciprocalOverlapTest(final int startA, final int endA, final int startB, final int endB, final double reciprocalOverlapEquals, final double reciprocalOverlapFalse) {
        final SVInterval intervalA = new SVInterval(0, startA, endA);
        final SVInterval intervalB1 = new SVInterval(0, startB, endB);
        Assert.assertEquals(SVIntervalUtils.reciprocalOverlap(intervalA, intervalB1), reciprocalOverlapEquals);
        Assert.assertTrue(SVIntervalUtils.hasReciprocalOverlap(intervalA, intervalB1, reciprocalOverlapEquals));
        Assert.assertTrue(SVIntervalUtils.hasReciprocalOverlap(intervalA, intervalB1, reciprocalOverlapEquals / 2));
        Assert.assertFalse(SVIntervalUtils.hasReciprocalOverlap(intervalA, intervalB1, reciprocalOverlapFalse));
        final SVInterval intervalB2 = new SVInterval(1, startB, endB);
        Assert.assertEquals(SVIntervalUtils.reciprocalOverlap(intervalA, intervalB2), 0.);
    }

    @Test(groups = "sv")
    public void convertIntervalsTest() {
        final SVInterval svInterval = new SVInterval(0, 100, 200);
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("seq0", 1000));
        final SimpleInterval simpleInterval = SVIntervalUtils.convertToSimpleInterval(svInterval, dictionary);
        Assert.assertEquals(simpleInterval, new SimpleInterval("seq0", 101, 200));
    }

}