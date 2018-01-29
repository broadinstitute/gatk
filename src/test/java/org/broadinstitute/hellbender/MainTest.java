package org.broadinstitute.hellbender;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import picard.cmdline.programgroups.Testing;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;

public final class MainTest extends CommandLineProgramTest {

    @Test(expectedExceptions = UserException.class)
    public void testCommandNotFoundThrows(){
        this.runCommandLine(new String[]{"Brain"});
    }

    @CommandLineProgramProperties(
            programGroup = Testing.class,
            summary = "OmitFromCommandLine test",
            oneLineSummary = "OmitFromCommandLine test",
            omitFromCommandLine = true)
    public static final class OmitFromCommandLineCLP extends CommandLineProgram {

        public static final int RETURN_VALUE = 1;

        @Override
        protected Object doWork() {
            return RETURN_VALUE;
        }
    }

    private static final class OmitFromCommandLineMain extends Main {
        @Override
        protected List<Class<? extends CommandLineProgram>> getClassList() {
            return Collections.singletonList(OmitFromCommandLineCLP.class);
        }
    }

    @Test
    public void testClpOmitFromCommandLine() {
        final OmitFromCommandLineMain main = new OmitFromCommandLineMain();
        final String clpName = "OmitFromCommandLineCLP";
        // test that the tool can be run from main correctly (returns non-null)
        Assert.assertEquals(main.instanceMain(new String[]{clpName}), OmitFromCommandLineCLP.RETURN_VALUE);
        // test that the usage is not shown if help is printed
        final String usage = captureStderr(() -> main.instanceMain(new String[]{"-h"}));
        Assert.assertFalse(usage.contains(clpName));
    }

}
