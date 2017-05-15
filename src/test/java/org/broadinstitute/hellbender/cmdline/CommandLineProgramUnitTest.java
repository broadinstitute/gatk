package org.broadinstitute.hellbender.cmdline;


import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.SpecialArgumentsCollection;
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
        assertContains(clp.getCommandLine(), SpecialArgumentsCollection.HELP_FULLNAME); //now it should be filled in
    }

    private static class ValidationFailer extends CommandLineProgram{
        public static final String ERROR1 = "first error";
        public static final String ERROR2 = "second error";

        @Override
        protected Object doWork() {
            return null;
        }

        @Override
        protected String[] customCommandLineValidation(){
            return new String[]{ERROR1, ERROR2};
        }
    }

    @Test
    public void testCustomValidationFailThrowsCommandLineException(){
        ValidationFailer clp = new ValidationFailer();
        try{
            clp.parseArgs(new String[] {});
            Assert.fail("Should have thrown an exception");
        } catch (final CommandLineException e){
            final String message = e.getMessage();
            assertContains(message, ValidationFailer.ERROR1);
            assertContains(message, ValidationFailer.ERROR2);
        }
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