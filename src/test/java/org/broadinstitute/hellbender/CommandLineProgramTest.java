package org.broadinstitute.hellbender;

import org.broadinstitute.hellbender.testutils.CommandLineProgramTester;

import java.io.File;
import java.util.List;

/**
 * Utility class for GATK CommandLine Program testing.
 */
public abstract class CommandLineProgramTest extends GATKBaseTest implements CommandLineProgramTester {

    /**
     * Returns the location of the resource directory. The default implementation points to the common directory for tools.
     */
    public static File getTestDataDir() {
        return new File("src/test/resources/org/broadinstitute/hellbender/tools/");
    }

    @Override
    public String getTestedToolName() {
        return getTestedClassName();
    }

    @Override
    public Object runCommandLine(final List<String> args) {
        return new Main().instanceMain(makeCommandLineArgs(args));
    }

    @Override
    public Object runCommandLine(final List<String> args, final String toolName) {
        return new Main().instanceMain(makeCommandLineArgs(args, toolName));
    }

}
