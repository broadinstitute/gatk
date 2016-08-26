package org.broadinstitute.hellbender;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.CommandLineProgramTestInterface;

import java.io.File;
import java.util.List;

/**
 * Utility class for GATK CommandLine Program testing.
 */
public abstract class CommandLineProgramTest extends BaseTest implements CommandLineProgramTestInterface {

    /**
     * Returns the location of the resource directory. The default implementation points to the common directory for tools.
     */
    public static File getTestDataDir() {
        return new File("src/test/resources/org/broadinstitute/hellbender/tools/");
    }

    public String getTestedToolName() {
        return getTestedClassName();
    }

    public Object runCommandLine(final List<String> args) {
        return new Main().instanceMain(makeCommandLineArgs(args));
    }

}
