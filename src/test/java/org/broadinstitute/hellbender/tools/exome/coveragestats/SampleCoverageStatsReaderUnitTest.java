package org.broadinstitute.hellbender.tools.exome.coveragestats;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Unit tests for {@link SampleCoverageStatsReader}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class SampleCoverageStatsReaderUnitTest extends BaseTest {

    private static final String TEST_SAMPLE_1 = "sample1";
    private static final String TEST_SAMPLE_2 = "sample2";

    @Test
    public void testRead() throws IOException {

        final File file = createTempFile("tcsru-test", ".tab");
        final PrintWriter writer = new PrintWriter(new FileWriter(file));
        writer.println(String.join("\t", SampleCoverageStats.SAMPLE_COLUMN_NAME,
                SampleCoverageStats.MEAN_COLUMN_NAME, SampleCoverageStats.VARIANCE_COLUMN_NAME));
        writer.println(String.join("\t", TEST_SAMPLE_1, "10", "100.19" ));

        writer.println(String.join("\t", TEST_SAMPLE_2, "101", "1010.19"));
        writer.close();
        final SampleCoverageStatsReader reader = new SampleCoverageStatsReader(file);
        final SampleCoverageStats stats1 = reader.readRecord();
        Assert.assertNotNull(stats1);
        Assert.assertEquals(stats1.mean, 10, 0.00000001);
        Assert.assertEquals(stats1.variance, 100.19, 0.0000001);
        Assert.assertEquals(stats1.sample, TEST_SAMPLE_1);
        final SampleCoverageStats stats2 = reader.readRecord();
        Assert.assertNotNull(stats2);
        Assert.assertEquals(stats2.mean, 101, 0.00000001);
        Assert.assertEquals(stats2.variance, 1010.19, 0.0000001);
        Assert.assertEquals(stats2.sample, TEST_SAMPLE_2);
        Assert.assertNull(reader.readRecord());
    }
}
