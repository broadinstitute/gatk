package org.broadinstitute.hellbender;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

public final class MainTest extends CommandLineProgramTest{

    @Test(expectedExceptions = UserException.class)
    public void testCommandNotFoundThrows(){
        this.runCommandLine(new String[]{"Brain"});
    }

    @Test
    public void testPrintReadsVersion() {
        String out = BaseTest.captureStderr(() -> Main.main(new String[]{"PrintReads", "--version"}));
        BaseTest.assertContains(out, "Version:");
    }
}
