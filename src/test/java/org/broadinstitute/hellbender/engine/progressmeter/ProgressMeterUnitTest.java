package org.broadinstitute.hellbender.engine.progressmeter;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class ProgressMeterUnitTest extends GATKBaseTest {

    @DataProvider(name = "UpdateIntervalTestData")
    public Object[][] getUpdateIntervalTestData() {
        return new Object[][] {
                // Seconds between logger updates, time waited, total number of records to process, and expected number of progress updates
                { 0.001, 5, 100},
                { 0.2, 5, 10000},
        };
    }

    @Test(dataProvider = "UpdateIntervalTestData")
    public void testUpdateInterval( final double secondsBetweenUpdates, final int waitTime, final long numRecords) throws InterruptedException {
        final ProgressMeter<Locatable> meter = new LocatableProgressMeter(secondsBetweenUpdates);
        //should take ~0 seconds...
        try {
            meter.start();
            Thread.sleep(10000);
            for (int i = 1; i <= numRecords; ++i) {
                meter.update(new SimpleInterval("1", 1, 1));
            }
            //now run the clock
            Thread.sleep(waitTime * 1000);
        } finally {
            meter.stop();
        }
        final double expectedUpdates = waitTime / secondsBetweenUpdates;
        final long updates = meter.numLoggerUpdates();
        final double allowableRatio = 1.0;
        BaseTest.assertEqualsDoubleSmart((updates - expectedUpdates) / expectedUpdates, expectedUpdates, allowableRatio,
                          "Wrong number of logger updates given secondsBetweenUpdates = " + secondsBetweenUpdates + " saw " + updates + " but expected between ");

        Assert.assertEquals(meter.getNumRecordsProcessed(), numRecords);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInvalidUpdateInterval() {
        final ProgressMeter<Locatable> meter = new LocatableProgressMeter(0.0);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testMustStartBeforeUpdate() {
        ProgressMeter<Locatable> pm = new LocatableProgressMeter();
        pm.update(new SimpleInterval("1",1,1));
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testCantStartTwice() {
        ProgressMeter<Locatable> pm = new LocatableProgressMeter();
        pm.start();
        pm.update(new SimpleInterval("1",1,1));
        pm.start();
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testCantStopBeforeStart() {
        ProgressMeter<Locatable> pm = new LocatableProgressMeter();
        pm.stop();
    }

    @Test
    public void testStartedAndStopped() {
        ProgressMeter<Locatable> pm = new LocatableProgressMeter();
        Assert.assertFalse(pm.started());
        Assert.assertFalse(pm.stopped());
        pm.start();
        Assert.assertTrue(pm.started());
        Assert.assertFalse(pm.stopped());
        pm.stop();
        Assert.assertTrue(pm.started());
        Assert.assertTrue(pm.stopped());
    }

}
