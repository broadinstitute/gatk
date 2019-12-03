package org.broadinstitute.hellbender.utils.variant;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;


public class ConvertToLocalAllelesIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testSomeNastyMultis(){
        ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(new File("multi-alleleics.04C29170.vcf"))
                .addOutput(new File("out.04C29170.vcf"));

        runCommandLine(args);
    }
}