package org.broadinstitute.hellbender.tools.walkers;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;


public class PrintBaseQualitiesIntegrationTest extends CommandLineProgramTest {
    @Test
    public void test(){
        final List<String> args = Arrays.asList(
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20");
        runCommandLine(args);
        final int d = 3;
    }



}