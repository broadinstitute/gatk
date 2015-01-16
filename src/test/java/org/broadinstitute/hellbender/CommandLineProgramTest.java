package org.broadinstitute.hellbender;

import org.broadinstitute.hellbender.utils.test.BaseTest;

import java.io.File;
import java.util.*;

/**
 * Utility class for CommandLine Program testing.
 */
public abstract class CommandLineProgramTest extends BaseTest {

    /**
     * Returns the location of the resource directory. The default implementation points to the common directory for tools.
     * Override if needed.
     */
    public static File getTestDataDir(){
        return new File("src/test/resources/org/broadinstitute/hellbender/tools/");
    }

    /**
     * Returns the name of the command line program name (ie the tool to use).
     * The default implementation takes the simple name of the test class and removes the trailing "Test".
     * Override if needed.
     */
    public String getCommandLineProgramName(){
        if (getClass().getSimpleName().contains("IntegrationTest"))
            return getClass().getSimpleName().replaceAll("IntegrationTest$", "");
        else
            return getClass().getSimpleName().replaceAll("Test$", "");
    }

    /**
     * For testing support.  Given a name of a Main CommandLineProgram and it's arguments, builds the arguments appropriate for calling the
     * program through Main
     *
     * @param args
     * @return String[] of command line arguments
     */
    public String[] makeCommandLineArgs(final List<String> args) {
        List<String> curatedArgs = injectDefaultVerbosity(args);
        final String[] commandLineArgs = new String[curatedArgs.size() + 1];
        commandLineArgs[0] = getCommandLineProgramName();
        int i = 1;
        for (final String arg : curatedArgs) {
            commandLineArgs[i++] = arg;
        }
        return commandLineArgs;
    }

    /**
     * Look for VERBOSITY argument; if not found, supply a default value that minimizes the amount of logging output.
     */
    private List<String> injectDefaultVerbosity(final List<String> args) {
        for (String arg : args) {
            if (arg.equalsIgnoreCase("--VERBOSITY")) return args;
        }
        List<String> argsWithVerbosity = new ArrayList<>(args);
        argsWithVerbosity.add("--VERBOSITY");
        argsWithVerbosity.add("ERROR");
        return argsWithVerbosity;
    }

    public String[] makeCommandLineArgs(final String[] args) {
        return makeCommandLineArgs(Arrays.asList(args));
    }

    public Object runCommandLine(final List<String> args) {
        return new Main().instanceMain(makeCommandLineArgs(args));
    }

    public Object runCommandLine(final String[] args) {
        return new Main().instanceMain(makeCommandLineArgs(args));
    }
}
