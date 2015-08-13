package org.broadinstitute.hellbender;

import com.google.common.base.Strings;
import org.broadinstitute.hellbender.dev.BunnyLog;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
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

    private static final String PROPERTY_DATAFLOW_RUNNER = "dataflowRunner";
    private static final String ENV_DATAFLOW_RUNNER = "DATAFLOW_RUNNER";

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

        // global toggle for BunnyLog output.
        BunnyLog.setEnabled(false);

        for (String arg : args) {
            if (arg.equalsIgnoreCase("--VERBOSITY")) return args;
        }
        List<String> argsWithVerbosity = new ArrayList<>(args);
        argsWithVerbosity.add("--VERBOSITY");
        argsWithVerbosity.add("INFO");
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

    /**
     * @return The unqualified class name of the dataflow runner if specified by either the
     * <code>dataflowRunner</code> system property or the <code>DATAFLOW_RUNNER</code> environment variable,
     * or <code>null</code> if no runner is specified.
     */
    public static String getExternallySpecifiedRunner() {
        String dataflowRunnerProperty = System.getProperty(PROPERTY_DATAFLOW_RUNNER);
        String dataflowRunnerEnv = System.getenv(ENV_DATAFLOW_RUNNER);
        if (!Strings.isNullOrEmpty(dataflowRunnerProperty)) {
            return dataflowRunnerProperty;
        } else if (!Strings.isNullOrEmpty(dataflowRunnerEnv)) {
            return dataflowRunnerEnv;
        } else {
            return null;
        }
    }

    /**
     * Adds arguments to the given <code>args</code> to set the dataflow
     * runner if specified by either the <code>dataflowRunner</code> system property
     * or the <code>DATAFLOW_RUNNER</code> environment variable (which
     * should be the unqualified class name of the runner).
     */
    public static void addDataflowRunnerArgs(ArgumentsBuilder args) {
        String dataflowRunner = getExternallySpecifiedRunner();
        if (!Strings.isNullOrEmpty(dataflowRunner)) {
            args.add("--runner");
            args.add(DataflowCommandLineProgram.getRunnerTypeName(dataflowRunner));
        }
    }

    /**
     * Adds arguments to the given <code>args</code> to set the dataflow
     * runner if specified by either the <code>dataflowRunner</code> system property
     * or the <code>DATAFLOW_RUNNER</code> environment variable (which
     * should be the unqualified class name of the runner).
     */
    public static void addDataflowRunnerArgs(List<String> args) {
        String dataflowRunner = getExternallySpecifiedRunner();
        if (!Strings.isNullOrEmpty(dataflowRunner)) {
            args.add("--runner");
            args.add(DataflowCommandLineProgram.getRunnerTypeName(dataflowRunner));
        }
    }
}
