package org.broadinstitute.hellbender;

import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.runtime.ProcessSettings;
import org.broadinstitute.hellbender.utils.test.CommandLineProgramTester;
import org.testng.Assert;

import java.io.File;
import java.util.Arrays;
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

}
