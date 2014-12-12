package org.broadinstitute.gatk;

import java.util.Arrays;
import java.util.List;

/**
 * Utility class for CommandLine Program testing.
 */
public abstract class CommandLineProgramTest {
    public abstract String getCommandLineProgramName();

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

    public int runCommandLine(final List<String> args) {
        return new Main().instanceMain(makeCommandLineArgs(args));
    }

    public int runCommandLine(final String[] args) {
        return new Main().instanceMain(makeCommandLineArgs(args));
    }
}
