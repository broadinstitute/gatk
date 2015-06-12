package org.broadinstitute.hellbender.cmdline;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.annotations.Test;

public class CommandLineProgramUnitTest {

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testExceptionInDoWorkNotSurpressedByExceptionInOnShutdown() {

        CommandLineProgram clp = new CommandLineProgram(){

            @Override
            protected Object doWork() {
                throw new UserException.BadInput("fail");
            }

            @Override
            protected void onShutdown(){
                throw new GATKException("another fail");
            }
        };
        clp.runTool();
    }

    @Test(expectedExceptions = GATKException.class)
    public void testExceptionInOnShutdownPropagatesIfNoEarlierExceptionOccurs(){
        CommandLineProgram clp = new CommandLineProgram(){

            @Override
            protected Object doWork() {
                return 0;
            }

            @Override
            protected void onShutdown(){
                throw new GATKException("fail");
            }
        };
        clp.runTool();
    }

}