package org.broadinstitute.hellbender;

import org.broadinstitute.hellbender.utils.test.BaseTest;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

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
     * @param args
     * @return String[] of command line arguments
     */
    public String[] makeCommandLineArgs(final List<String> args) {
        List<String> curatedArgs = injectDefaultVerbosity(args);
        final String[] commandLineArgs = new String[curatedArgs.size() + 1];
        commandLineArgs[0] = getTestedClassName();
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

    /**
     * This finds the first file that is in the same directory as the given prefix file, which is not the given file,
     * and which starts with the same filename as the given file.
     * @throws IOException
     */
    public static File findDataflowOutput(File pathPrefix) throws IOException{
         Optional<File> outputFile = Files.list(pathPrefix.toPath().getParent())
                .filter(f -> f.getFileName().toString().startsWith(pathPrefix.getName()))
                .filter(f -> !f.equals(pathPrefix.toPath()))
                .map(Path::toFile)
                .findFirst();

        File file =  outputFile.orElseThrow(() -> new IOException("No dataflow output file find for prefix: "+ pathPrefix.getAbsolutePath()));
        file.deleteOnExit();
        return file;
    }
}
