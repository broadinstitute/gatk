package org.broadinstitute.hellbender.utils.help;

import org.broadinstitute.barclay.argparser.ClassFinder;
import org.broadinstitute.barclay.help.BashTabCompletionDoclet;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.testng.annotations.Test;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Smoke test to run tab completion generation on a subset of classes to make sure it doesn't regress.
 */
public class TabCompletionIntegrationTest extends CommandLineProgramTest {

    /**
     * Entry point for manually running the gatkTabComplete process on a subset of packages from within GATK.
     */
    final private static List<String> tabCompletionTestPackages;

    // Static block to initialize tabCompletionTestPackages:
    static {
        final ClassFinder classFinder = new ClassFinder();
        classFinder.find("org.broadinstitute.hellbender", CommandLineProgram.class);
        tabCompletionTestPackages = Collections.unmodifiableList(
                classFinder.getClasses().stream().map(Class::getName).collect(Collectors.toList())
        );
    }

    @Test
    public static void tabCompleteSmokeTest() throws IOException, InterruptedException {
        final File tabCompletionTestTarget = createTempDir("tabCompletionTest");

        // Setup rote input arguments:
        final List<String> argList = Arrays.asList(

                "javadoc",

                "-doclet", BashTabCompletionDoclet.class.getName(),

                "-docletpath", System.getProperty("java.class.path"),
                "-sourcepath", "src/main/java",
                "-d", tabCompletionTestTarget.getAbsolutePath(), // directory must exist
                "-cp", System.getProperty("java.class.path"),

                "-use-default-templates",

                "-output-file-extension", "sh",
                "-index-file-extension", "sh",
                "-absolute-version", "0.0-001",
                "-build-timestamp", new SimpleDateFormat("dd-mm-yyyy hh:mm:ss").format( new Date() ),

                "-caller-script-name", "gatk",

                "-caller-pre-legal-args", "--help --list --dry-run --java-options",
                "-caller-pre-arg-val-types", "null null null String",
                "-caller-pre-mutex-args", "--help;list,dry-run,java-options --list;help,dry-run,java-options",
                "-caller-pre-alias-args", "--help;-h",
                "-caller-pre-arg-min-occurs", "0 0 0 0",
                "-caller-pre-arg-max-occurs", "1 1 1 1",

                "-caller-post-legal-args", "--spark-runner --spark-master --cluster --dry-run --java-options --conf --driver-memory --driver-cores --executor-memory --executor-cores --num-executors",
                "-caller-post-arg-val-types", "String String String null String file int int int int int",
                "-caller-post-mutex-args", "",
                "-caller-post-alias-args", "",
                "-caller-post-arg-min-occurs", "0 0 0 0 0 0 0 0 0 0",
                "-caller-post-arg-max-occurs", "1 1 1 1 1 1 1 1 1 1",

                "-verbose"
        );

        // Point the javadoc at our packages:
        final List<String> docArgList = new ArrayList<>();
        docArgList.addAll(argList);
        docArgList.addAll(tabCompletionTestPackages);

        // This is  smoke test; we just want to make sure it doesn't blow up

        // Run this as a process, not through Java itself:
        final ProcessController processController = new ProcessController();
        runProcess(processController, docArgList.toArray(new String[] {}) );
    }

}
