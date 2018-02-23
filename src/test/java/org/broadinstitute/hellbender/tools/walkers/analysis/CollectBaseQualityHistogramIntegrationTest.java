package org.broadinstitute.hellbender.tools.walkers.analysis;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

import java.io.File;

import static org.testng.Assert.*;

/**
 * Created by tsato on 2/23/18.
 */
public class CollectBaseQualityHistogramIntegrationTest extends CommandLineProgramTest {

    @Test
    public void test(){
        final File outputIntervalFile = createTempFile("temp", ".bed");
        final String[] args = {
                "-I", NA12878_chr17_1k_BAM,
                "-R", v37_chr17_1Mb_Reference,
                "-O", outputIntervalFile.getAbsolutePath()
        };

        runCommandLine(args);
    }
}