package org.broadinstitute.hellbender.engine.dataflow;

import com.google.common.base.Strings;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

public class DataflowCommandLineProgramTest extends CommandLineProgramTest {

    private static final String PROPERTY_DATAFLOW_RUNNER = "dataflowRunner";
    private static final String ENV_DATAFLOW_RUNNER = "DATAFLOW_RUNNER";

    /**
     * @return an arguments List containing the apiKey, project, and staging arguments
     *         populated from environment variables as defined in {@link #getGCPTestApiKey},
     *         {@link #getGCPTestProject}, and {@link #getGCPTestStaging}, suitable
     *         for use in a hellbender command line.
     */
    public static List<String> getStandardGCPArgumentsFromEnvironment() {
        return Arrays.asList("--apiKey", getGCPTestApiKey(),
                             "--project", getGCPTestProject(),
                             "--staging", getGCPTestStaging());
    }

    /**
     * This finds the first file that is in the same directory as the given prefix file, which is not the given file,
     * and which starts with the same filename as the given file.
     * @throws java.io.IOException
     */
    public static File findDataflowOutput(File pathPrefix) throws IOException {
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
