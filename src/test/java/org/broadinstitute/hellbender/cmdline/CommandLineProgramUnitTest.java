package org.broadinstitute.hellbender.cmdline;


import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.junit.Assert;
import org.testng.annotations.Test;

public class CommandLineProgramUnitTest extends BaseTest {

    @Test
    public void testGetUsage(){
        final CommandLineProgram clp = getClp();
        String usage = clp.getUsage();
        BaseTest.assertContains(usage, "Usage:");
    }

    @Test
    public void testGetCommandLine(){
        final CommandLineProgram clp = getClp();
        Assert.assertNull(clp.getCommandLine()); //should be null since no args were specified
        Assert.assertFalse(clp.parseArgs(new String[]{"--" + SpecialArgumentsCollection.HELP_FULLNAME}));
        assertContains(clp.getCommandLine(), SpecialArgumentsCollection.HELP_FULLNAME); //now it should be filed in
    }

    private static CommandLineProgram getClp() {
        return new CommandLineProgram() {
            @Override
            protected Object doWork() {
                return null;
            }
        };
    }
}