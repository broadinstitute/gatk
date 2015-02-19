package org.broadinstitute.hellbender;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

public class MainTest extends CommandLineProgramTest{
    @Test
    public void test(){
        //this just tests that the program does not crash
        this.runCommandLine(new String[]{"Brain"});
    }

    @Test
    public void testPrintReadsVersion() {
        String out = BaseTest.captureStderr(() -> Main.main(new String[]{"PrintReads", "--version"}));
        BaseTest.assertContains(out, "Version:");
    }

    @Test
    public void testUserException() {
        String out = BaseTest.captureStderr(() -> Main.main(new String[]{"PrintReads", "blargle"}));
        BaseTest.assertContains(out, "USER ERROR");
    }

}
