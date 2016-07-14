package org.broadinstitute.hellbender.tools.exome.coveragestats;

import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Unit tests for {@link SampleCoverageStatsReader}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class TargetCoverageStatsWriterUnitTest extends BaseTest {

    private static final Target TEST_TARGET_1 = new Target("target1", new SimpleInterval("chr1", 1, 100));
    private static final Target TEST_TARGET_2 = new Target("target2", new SimpleInterval("chr2", 300, 400));

    @Test
    public void testWithIntervals() throws IOException {
        final File file = createTempFile("tcsw-test", ".tab");
        final TargetCoverageStatsWriter writer = new TargetCoverageStatsWriter(file, true);
        writer.writeRecord(new TargetCoverageStats(TEST_TARGET_1, 10, 100, 10));
        writer.writeRecord(new TargetCoverageStats(TEST_TARGET_2, 101, 1212.12, 30));
        writer.close();
        final TargetCoverageStatsReader reader = new TargetCoverageStatsReader(file);
        final TargetCoverageStats stats1 = reader.readRecord();
        Assert.assertNotNull(stats1);
        Assert.assertEquals(stats1.target, TEST_TARGET_1);
        Assert.assertEquals(stats1.target.getInterval(), TEST_TARGET_1.getInterval());
        Assert.assertEquals(stats1.mean, 10, 0.0000001);
        Assert.assertEquals(stats1.variance, 100, 0.0000001);
        Assert.assertEquals(stats1.interquartileRange, 10, 0.0000001);
        final TargetCoverageStats stats2 = reader.readRecord();
        Assert.assertNotNull(stats2);
        Assert.assertEquals(stats2.target, TEST_TARGET_2);
        Assert.assertEquals(stats2.target.getInterval(), TEST_TARGET_2.getInterval());
        Assert.assertEquals(stats2.mean, 101, 0.0000001);
        Assert.assertEquals(stats2.variance, 1212.12, 0.0000001);
        Assert.assertEquals(stats2.interquartileRange, 30, 0.0000001);
        Assert.assertNull(reader.readRecord());
    }

    @Test
    public void testWithoutIntervals() throws IOException {
        final File file = createTempFile("tcsw-test", ".tab");
        final TargetCoverageStatsWriter writer = new TargetCoverageStatsWriter(file, false);
        writer.writeRecord(new TargetCoverageStats(TEST_TARGET_1, 10, 100, 10));
        writer.writeRecord(new TargetCoverageStats(TEST_TARGET_2, 101, 1212.12, 30));
        writer.close();
        final TargetCoverageStatsReader reader = new TargetCoverageStatsReader(file);
        final TargetCoverageStats stats1 = reader.readRecord();
        Assert.assertNotNull(stats1);
        Assert.assertEquals(stats1.target, TEST_TARGET_1);
        Assert.assertNull(stats1.target.getInterval());
        Assert.assertEquals(stats1.mean, 10, 0.0000001);
        Assert.assertEquals(stats1.variance, 100, 0.0000001);
        Assert.assertEquals(stats1.interquartileRange, 10, 0.0000001);
        final TargetCoverageStats stats2 = reader.readRecord();
        Assert.assertNotNull(stats2);
        Assert.assertEquals(stats2.target, TEST_TARGET_2);
        Assert.assertNull(stats2.target.getInterval());
        Assert.assertEquals(stats2.mean, 101, 0.0000001);
        Assert.assertEquals(stats2.variance, 1212.12, 0.0000001);
        Assert.assertEquals(stats2.interquartileRange, 30, 0.0000001);
        Assert.assertNull(reader.readRecord());
    }
}
