package org.broadinstitute.hellbender.tools.walkers.misc;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;

import static org.testng.Assert.*;

/**
 * Created by tsato on 1/24/18.
 */
public class CreateHomopolymerIntervalsIntegrationTest extends CommandLineProgramTest{
    final File outputIntervalFile = createTempFile("temp", ".intervals");

    @Test
    public void test(){
        final String[] args = {
                "-R", hg19MiniReference,
                "-O", outputIntervalFile.getAbsolutePath()
        };

        runCommandLine(args);
        int d = 3;
    }


}