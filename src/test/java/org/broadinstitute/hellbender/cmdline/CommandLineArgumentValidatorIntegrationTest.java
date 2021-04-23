package org.broadinstitute.hellbender.cmdline;

import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineArgumentValidatorMain;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.Test;

public class CommandLineArgumentValidatorIntegrationTest extends GATKBaseTest {

    @Test
    public void testArgumentValidatorPositive() {
        CommandLineArgumentValidatorMain.main(new String[] {
                "PrintReads",
                "-I",
                "filesDontNeedToExistForCommandLineValidation.bam",
                "-O",
                "filesDontNeedToExistForCommandLineValidation.bam"
        });
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testArgumentValidatorNegative() {
        // force a command line parsing exception
        CommandLineArgumentValidatorMain.main(new String[] {
                "PrintReads",
                "-unrecognizedOption"
        });
    }

}
