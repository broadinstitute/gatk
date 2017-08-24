package org.broadinstitute.hellbender;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.junit.Assert;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;

public final class MainTest extends CommandLineProgramTest {

    @Test(expectedExceptions = UserException.class)
    public void testCommandNotFoundThrows(){
        this.runCommandLine(new String[]{"Brain"});
    }


    /** Main class for testing a single CLP */
    private static final class SingleClpMain extends Main {

        private Class<? extends CommandLineProgram> clazz;

        private SingleClpMain(final Class<? extends CommandLineProgram> clazz) {
            this.clazz = clazz;
        }

        @Override
        protected List<String> getPackageList() {
            return Collections.emptyList();
        }

        @Override
        protected List<Class<? extends CommandLineProgram>> getClassList() {
            return Collections.singletonList(clazz);
        }
    }

    @CommandLineProgramProperties(
            programGroup = TestProgramGroup.class,
            summary = "OmitFromCommanLine test",
            oneLineSummary = "OmitFromCommanLine test",
            omitFromCommandLine = true)
    public static final class OmitFromCommanLineCLP extends CommandLineProgram {

        @Override
        protected Object doWork() {
            return 1;
        }
    }

    @Test
    public void testClpOmitFromCommandLine() {
        final SingleClpMain main = new SingleClpMain(OmitFromCommanLineCLP.class);
        // test that the tool can be run from main correctly (returns non-null)
        Assert.assertNotNull(main.instanceMain(new String[]{"OmitFromCommanLineCLP"}));
        // test that the usage is not shown if help is printed
        final String usage = captureStderr(() -> main.instanceMain(new String[]{"-h"}));
        Assert.assertFalse(usage.contains("OmitFromCommanLineCLP"));
    }

}
