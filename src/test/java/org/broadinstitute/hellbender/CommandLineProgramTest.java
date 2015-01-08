package org.broadinstitute.hellbender;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * Utility class for CommandLine Program testing.
 */
public abstract class CommandLineProgramTest {

    /**
     * Returns the location of the resource directory. The default implementation points to the common directory for tools.
     * Override if needed.
     */
    protected File getTestDataDir(){
        return new File("src/test/resources/org/broadinstitute/hellbender/tools/");
    }

    /**
     * Returns the name of the command line program name (ie the tool to use).
     * The default implementation takes the simple name of the test class and removes the trailing "Test".
     * Override if needed.
     */
    public String getCommandLineProgramName(){
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
        final String[] commandLineArgs = new String[args.size() + 1];
        commandLineArgs[0] = getCommandLineProgramName();
        int i = 1;
        for (final String arg : args) {
            commandLineArgs[i++] = arg;
        }
        return commandLineArgs;
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
