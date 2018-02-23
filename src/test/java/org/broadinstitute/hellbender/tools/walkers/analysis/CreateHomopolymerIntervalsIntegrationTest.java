package org.broadinstitute.hellbender.tools.walkers.analysis;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Created by tsato on 1/24/18.
 */
public class CreateHomopolymerIntervalsIntegrationTest extends CommandLineProgramTest{
    @Test
    public void test() {
        final File outputIntervalFile = createTempFile("temp", ".bed");
        final String[] args = {
                "-R", hg19MiniReference,
                "-O", outputIntervalFile.getAbsolutePath()
        };

        runCommandLine(args);
        int d = 3;
    }

    @Test
    public void testEndOfChromosome17() {
        final File outputIntervalFile = createTempFile("temp", ".bed");
        final String[] args = {
                "-R", hg19EndOfChromosome17,
                "-O", outputIntervalFile.getAbsolutePath(),
        };

        runCommandLine(args);
        Assert.assertEquals(outputIntervalFile.length(), 0);
    }



}