package org.broadinstitute.hellbender;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.security.Permission;
import java.util.Collections;
import java.util.List;

public final class MainTest extends CommandLineProgramTest {

    @Test(expectedExceptions = UserException.class)
    public void testCommandNotFoundThrows(){
        this.runCommandLine(new String[]{"Brain"});
    }

    @CommandLineProgramProperties(
            programGroup = TestProgramGroup.class,
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


    private static final class ExitNotAllowedException extends SecurityException {
        private static final long serialVersionUID = 1L;
        final int status;

        ExitNotAllowedException(int status) {
            this.status = status;
        }
    }

    private static final class ThrowOnExitSecurityManager extends SecurityManager {

        @Override
        public void checkPermission(Permission perm) {
            // allow anything.
        }

        @Override
        public void checkPermission(Permission perm, Object context) {
            // allow anything.
        }

        @Override
        public void checkExit(int status) {
            super.checkExit(status);
            // always throw
            throw new ExitNotAllowedException(status);
        }
    }

    @Test(singleThreaded = true)
    public void testMainErrorWithoutStackTrace() {
        final SecurityManager backup = System.getSecurityManager();
        try {
            System.setSecurityManager(new ThrowOnExitSecurityManager());
            new Main().mainEntry(new String[]{"PrintReadsW"});
            Assert.fail("Should never reach here");
        } catch (ExitNotAllowedException e) {
            // does exist as if it is an user exception
            Assert.assertEquals(e.status, Main.USER_EXCEPTION_EXIT_VALUE);
        } finally {
            System.setSecurityManager(backup);
        }
    }
    @Test(singleThreaded = true)
    public void testNonZeroPicardReturnValue() {
        final SecurityManager backup = System.getSecurityManager();
        try {
            System.setSecurityManager(new ThrowOnExitSecurityManager());
            new Main().mainEntry(new String[]{"ExtractSequences"});
            Assert.fail("Should never reach here");
        } catch (final ExitNotAllowedException e) {
            Assert.assertEquals(e.status, Main.PICARD_TOOL_EXCEPTION);
        } finally {
            System.setSecurityManager(backup);
        }
    }

}
