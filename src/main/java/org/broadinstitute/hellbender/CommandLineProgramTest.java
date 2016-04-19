package org.broadinstitute.hellbender;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.logging.BunnyLog;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Utility class for CommandLine Program testing.
 */
public abstract class CommandLineProgramTest extends BaseTest {

    /**
     * Returns the location of the resource directory. The default implementation points to the common directory for tools.
     */
    public static File getTestDataDir(){
        return new File("src/test/resources/org/broadinstitute/hellbender/tools/");
    }


    /**
     * For testing support.  Given a name of a Main CommandLineProgram and it's arguments, builds the arguments appropriate for calling the
     * program through Main
     *
     * @param args List<String> of command line arguments
     * @return String[] of command line arguments
     */
    public String[] makeCommandLineArgs(final List<String> args) {
        return makeCommandLineArgs(args, getTestedClassName());
    }

    /**
     * For testing support.  Given a name of a Main CommandLineProgram and it's arguments, builds the arguments appropriate for calling the
     * program through Main
     *
     * @param args List<String> of command line arguments
     * @param toolname name of the tool to test
     * @return String[] of command line arguments
     */
    public String[] makeCommandLineArgs(final List<String> args, final String toolname) {
        List<String> curatedArgs = injectDefaultVerbosity(args);
        final String[] commandLineArgs = new String[curatedArgs.size() + 1];
        commandLineArgs[0] = toolname;
        int i = 1;
        for (final String arg : curatedArgs) {
            commandLineArgs[i++] = arg;
        }
        return commandLineArgs;
    }

    /**
     * Look for --verbosity argument; if not found, supply a default value that minimizes the amount of logging output.
     */
    private List<String> injectDefaultVerbosity(final List<String> args) {

        // global toggle for BunnyLog output.
        BunnyLog.setEnabled(false);

        for (String arg : args) {
            if (arg.equalsIgnoreCase("--" + StandardArgumentDefinitions.VERBOSITY_NAME) || arg.equalsIgnoreCase("-" + StandardArgumentDefinitions.VERBOSITY_NAME)) {
                return args;
            }
        }
        List<String> argsWithVerbosity = new ArrayList<>(args);
        argsWithVerbosity.add("--" + StandardArgumentDefinitions.VERBOSITY_NAME);
        argsWithVerbosity.add(Log.LogLevel.ERROR.name());
        return argsWithVerbosity;
    }

    public Object runCommandLine(final List<String> args) {
        return new Main().instanceMain(makeCommandLineArgs(args));
    }

    public Object runCommandLine(final String[] args) {
        return runCommandLine(Arrays.asList(args));
    }

    public Object runCommandLine(final ArgumentsBuilder args) {
        return runCommandLine(args.getArgsList());
    }

}
