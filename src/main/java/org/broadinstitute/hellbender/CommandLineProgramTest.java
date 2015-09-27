package org.broadinstitute.hellbender;

import com.google.common.base.Strings;
import com.google.common.io.CharStreams;
import org.apache.commons.lang.StringUtils;
import org.apache.spark.launcher.SparkLauncher;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.logging.BunnyLog;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
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
     * @param args List<String> of command line arguments
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

    public String buildSparkJarIfMissing() throws IOException, InterruptedException {
        // Check for the sparkJar, if missing make it
        Runtime r = Runtime.getRuntime();
        // Use gradle properties to find the jar name
        // TODO: we could try to use git describe --always
        Process jarP = r.exec("gradle properties");
        CharStreams.toString(new InputStreamReader(jarP.getErrorStream()));
        final String[] properties = CharStreams.toString(new InputStreamReader(jarP.getInputStream())).split("\\r?\\n");
        String version = "";
        for (String s : properties) {
            int idx = s.indexOf("version: ");
            if (idx != -1) {
                version = s.substring(9);
            }
        }

        if (version.isEmpty()) {
            throw new GATKException("unable to find version from \"gradle properties\"");
        }

        System.out.println("**** " + version);
        String fullJarName = "hellbender_branch-all-" + version + "-spark.jar";
        String jarPath = "build/libs/" + fullJarName;

        Process jarSearchP = r.exec("test -f " + jarPath);
        //CharStreams.toString(new InputStreamReader(jarSearchP.getInputStream())); // Drain the streams.
        //CharStreams.toString(new InputStreamReader(jarSearchP.getErrorStream())); // Drain the streams.
        if (jarSearchP.waitFor() != 0) {
            System.out.println("jar not found");
            // Build the jar
            Process shadowJarP = r.exec("gradle sparkJar");
            //CharStreams.toString(new InputStreamReader(shadowJarP.getInputStream())); // Drain the streams.
            //CharStreams.toString(new InputStreamReader(shadowJarP.getErrorStream())); // Drain the streams.
            int exitCode = shadowJarP.waitFor();
            if (exitCode != 0) {
                throw new GATKException("unable to build the Spark jar");
            }
        }
        return jarPath;
    }

    public int runLocalSpark(final List<String> args) throws IOException, InterruptedException {
        args.add("--sparkMaster");
        args.add(getSparkMasterURL());
        String path = new File("build/classes/main").getCanonicalPath();
        String jarLibPath = new File("build/lib/*").getCanonicalPath();
        String[] cmdArgs = makeCommandLineArgs(args);


        Process spark = new SparkLauncher()
                .setAppResource("")
                .setMainClass("org.broadinstitute.hellbender.Main")
                .setMaster(getSparkMasterURL()) //"spark://davidada-macbookpro2.roam.corp.google.com:7077")
                .setConf(SparkLauncher.EXECUTOR_EXTRA_CLASSPATH, path + ":" + jarLibPath)
                .setConf(SparkLauncher.DRIVER_EXTRA_CLASSPATH, path + ":" + jarLibPath)
                        .addAppArgs(cmdArgs)
                .launch();
        String stdout = CharStreams.toString(new InputStreamReader(spark.getInputStream()));// Drain the streams.
        String stderr = CharStreams.toString(new InputStreamReader(spark.getErrorStream()));// Drain the streams.
        int returnCode = spark.waitFor();
        System.out.println(stdout);
        System.err.println(stderr);
        return returnCode;
    }

    public int runSparkGoogleCommandLine(final List<String> args) throws IOException, InterruptedException {
        args.add("--sparkMaster");
        args.add("yarn-client");
        String[] cmdArgs = makeCommandLineArgs(args);

        String jarPath = buildSparkJarIfMissing();
        String dataProcCommand = "gcloud alpha dataproc jobs submit spark --cluster "
                + getGoogleSparkCluster() + " --class org.broadinstitute.hellbender.Main --jar "
                + jarPath + " " + StringUtils.join(cmdArgs, " ");

        Runtime r = Runtime.getRuntime();

        Process dataProc = r.exec(dataProcCommand);
        String stdout = CharStreams.toString(new InputStreamReader(dataProc.getInputStream()));
        String stderr = CharStreams.toString(new InputStreamReader(dataProc.getErrorStream()));
        int returnCode = dataProc.waitFor();
        System.out.println(stdout);
        System.out.println(stderr);
        return returnCode;
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
