package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import static org.testng.Assert.*;

public class MergeMutectStatsIntegrationTest extends CommandLineProgramTest {

    @Test
    public void simpleTest() {
        final File statsFile1 = createTempFile("stats1", ".stats");
        final File statsFile2 = createTempFile("stats2", ".stats");
        final File merged = createTempFile("merged", ".stats");

        final List<MutectStats> stats1 = Arrays.asList(new MutectStats(Mutect2Engine.CALLABLE_SITES_NAME, 20));
        final List<MutectStats> stats2 = Arrays.asList(new MutectStats(Mutect2Engine.CALLABLE_SITES_NAME, 30));

        MutectStats.writeToFile(stats1, statsFile1);
        MutectStats.writeToFile(stats2, statsFile2);

        final List<String> args = Arrays.asList("--" + Mutect2.MUTECT_STATS_SHORT_NAME, statsFile1.getAbsolutePath(),
                "--" + Mutect2.MUTECT_STATS_SHORT_NAME, statsFile2.getAbsolutePath(),
                "-O", merged.getAbsolutePath());
        runCommandLine(args);

        final List<MutectStats> mergedStats = MutectStats.readFromFile(merged);

        Assert.assertEquals(mergedStats.size(), 1);
        Assert.assertEquals(mergedStats.get(0).getStatistic(), Mutect2Engine.CALLABLE_SITES_NAME);
        Assert.assertEquals(mergedStats.get(0).getValue(), 50, 1e-9);
    }



}