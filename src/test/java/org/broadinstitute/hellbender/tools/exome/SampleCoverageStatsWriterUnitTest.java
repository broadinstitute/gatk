package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Unit tests for {@link SampleCoverageStatsWriter}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class SampleCoverageStatsWriterUnitTest extends BaseTest {

    private static final String TEST_SAMPLE_1 = "sample1";
    private static final String TEST_SAMPLE_2 = "sample2";

    @Test
    public void testWrite() throws IOException {
        final File file = createTempFile("scsw-test", ".tab");
        final SampleCoverageStatsWriter writer = new SampleCoverageStatsWriter(file);
        writer.writeRecord(new SampleCoverageStats(TEST_SAMPLE_1, 10, 101.1));
        writer.writeRecord(new SampleCoverageStats(TEST_SAMPLE_2, 100, 101.2));
        writer.close();
        final SampleCoverageStatsReader reader = new SampleCoverageStatsReader(file);
        final SampleCoverageStats stats1 = reader.readRecord();
        Assert.assertNotNull(stats1);
        Assert.assertEquals(stats1.mean, 10.0, 0.000000001);
        Assert.assertEquals(stats1.variance, 101.1, 0.000000001);
        final SampleCoverageStats stats2 = reader.readRecord();
        Assert.assertNotNull(stats2);
        Assert.assertEquals(stats2.mean, 100, 0.000000001);
        Assert.assertEquals(stats2.variance, 101.2, 0.000000001);
        Assert.assertNull(reader.readRecord());
    }
}
