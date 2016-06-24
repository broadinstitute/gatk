package org.broadinstitute.hellbender.tools.exome.coveragestats;

import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetTableColumn;
import org.broadinstitute.hellbender.tools.exome.coveragestats.TargetCoverageStats;
import org.broadinstitute.hellbender.tools.exome.coveragestats.TargetCoverageStatsReader;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Unit tests for {@link TargetCoverageStatsReader}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class TargetCoverageStatsReaderUnitTest extends BaseTest {

    private static final Target TEST_TARGET_1 = new Target("target_1", new SimpleInterval("chr1", 1, 100));
    private static final Target TEST_TARGET_2 = new Target("target_2", new SimpleInterval("chr2", 200, 300));

    @Test
    public void testReadWithIntervals() throws IOException {

        final File file = createTempFile("tcsru-test", ".tab");
        final PrintWriter writer = new PrintWriter(new FileWriter(file));
        writer.println(String.join("\t", TargetTableColumn.NAME.toString(), TargetTableColumn.CONTIG.toString(),
                TargetTableColumn.START.toString(), TargetTableColumn.END.toString(),
                TargetCoverageStats.MEAN_COLUMN_NAME, TargetCoverageStats.VARIANCE_COLUMN_NAME,
                TargetCoverageStats.INTERQUARTILE_RANGE_COLUMN_NAME));
        writer.println(String.join("\t", TEST_TARGET_1.getName(), TEST_TARGET_1.getInterval().getContig(),
                "" + TEST_TARGET_1.getStart(), "" + TEST_TARGET_1.getEnd(), "10", "100.19", "70.3" ));

        writer.println(String.join("\t", TEST_TARGET_2.getName(), TEST_TARGET_2.getInterval().getContig(),
                "" + TEST_TARGET_2.getStart(), "" + TEST_TARGET_2.getEnd(), "101", "1010.19", "703.23"));
        writer.close();
        final TargetCoverageStatsReader reader = new TargetCoverageStatsReader(file);
        final TargetCoverageStats stats1 = reader.readRecord();
        Assert.assertNotNull(stats1);
        Assert.assertEquals(stats1.mean, 10, 0.00000001);
        Assert.assertEquals(stats1.variance, 100.19, 0.0000001);
        Assert.assertEquals(stats1.interquartileRange, 70.3, 0.0000001);
        Assert.assertEquals(stats1.target.getName(), TEST_TARGET_1.getName());
        Assert.assertEquals(stats1.target.getInterval(), TEST_TARGET_1.getInterval());
        final TargetCoverageStats stats2 = reader.readRecord();
        Assert.assertNotNull(stats2);
        Assert.assertEquals(stats2.mean, 101, 0.00000001);
        Assert.assertEquals(stats2.variance, 1010.19, 0.0000001);
        Assert.assertEquals(stats2.interquartileRange, 703.23, 0.0000001);
        Assert.assertEquals(stats2.target.getName(), TEST_TARGET_2.getName());
        Assert.assertEquals(stats2.target.getInterval(), TEST_TARGET_2.getInterval());
        Assert.assertNull(reader.readRecord());
    }

    @Test
    public void testWithoutIntervals() throws IOException {

        final File file = createTempFile("tcsru-test", ".tab");
        final PrintWriter writer = new PrintWriter(new FileWriter(file));
        writer.println(String.join("\t", TargetTableColumn.NAME.toString(),
                TargetCoverageStats.MEAN_COLUMN_NAME, TargetCoverageStats.VARIANCE_COLUMN_NAME,
                TargetCoverageStats.INTERQUARTILE_RANGE_COLUMN_NAME));
        writer.println(String.join("\t", TEST_TARGET_1.getName(), "10", "100.19", "70.3" ));

        writer.println(String.join("\t", TEST_TARGET_2.getName(), "101", "1010.19", "703.23"));
        writer.close();
        final TargetCoverageStatsReader reader = new TargetCoverageStatsReader(file);
        final TargetCoverageStats stats1 = reader.readRecord();
        Assert.assertNotNull(stats1);
        Assert.assertEquals(stats1.mean, 10, 0.00000001);
        Assert.assertEquals(stats1.variance, 100.19, 0.0000001);
        Assert.assertEquals(stats1.interquartileRange, 70.3, 0.0000001);
        Assert.assertEquals(stats1.target.getName(), TEST_TARGET_1.getName());
        Assert.assertNull(stats1.target.getInterval());
        final TargetCoverageStats stats2 = reader.readRecord();
        Assert.assertNotNull(stats2);
        Assert.assertEquals(stats2.mean, 101, 0.00000001);
        Assert.assertEquals(stats2.variance, 1010.19, 0.0000001);
        Assert.assertEquals(stats2.interquartileRange, 703.23, 0.0000001);
        Assert.assertEquals(stats2.target.getName(), TEST_TARGET_2.getName());
        Assert.assertNull(stats2.target.getInterval());
        Assert.assertNull(reader.readRecord());
    }

}
