package org.broadinstitute.hellbender.tools.picard.vcf;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by lichtens on 7/11/17.
 */
public class VcfToIntervalListIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_RESOURCE_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome/orientationbiasvariantfilter/");
    public static final String smallM2VcfMore = TEST_RESOURCE_DIR.getAbsolutePath() + "/small_m2_more_variants.vcf";

    @Test
    public void testExcludingFiltered() throws IOException {
        final File outputFile = File.createTempFile("vcftointervallist_", ".interval_list");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(smallM2VcfMore);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final List<Interval> intervals = IntervalList.fromFile(outputFile).getIntervals();

        // 11 total, 10 defined unique intervals (two variants are adjacent), one is filtered in INFO and
        // one is filtered in FORMAT FT, but only INFO counts
        Assert.assertEquals(intervals.size(), 11 - 1 - 1);
    }

    @Test
    public void testIncludingFiltered() throws IOException {
        final File outputFile = File.createTempFile("vcftointervallist_incl_", ".interval_list");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        arguments.add(smallM2VcfMore);
        arguments.add("-" + VcfToIntervalList.INCLUDE_FILTERED_SHORT_NAME);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final List<Interval> intervals = IntervalList.fromFile(outputFile).getIntervals();

        // 11 total, 10 defined unique intervals (two variants are adjacent)
        Assert.assertEquals(intervals.size(), 11 - 1 );
    }
}
