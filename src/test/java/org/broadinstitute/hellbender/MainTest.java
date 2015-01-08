package org.broadinstitute.hellbender;

import org.testng.annotations.Test;

public class MainTest extends CommandLineProgramTest{

    @Test
    public void test(){
        //this just tests that the program does not crash
        this.runCommandLine(new String[]{"Brain"});
    }
}
