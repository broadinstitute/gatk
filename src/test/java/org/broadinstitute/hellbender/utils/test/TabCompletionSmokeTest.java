package org.broadinstitute.hellbender.utils.test;

import htsjdk.samtools.util.ProcessExecutor;
import org.broadinstitute.barclay.help.BashTabCompletionDoclet;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ClassFinder;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.help.GATKHelpDoclet;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.broadinstitute.hellbender.utils.runtime.ProcessSettings;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Smoke test to run tab completion generation on a subset of classes to make sure it doesn't regress.
 */
public class TabCompletionSmokeTest extends CommandLineProgramTest {

    /**
     * Entry point for manually running the gatkTabComplete process on a subset of packages from within GATK.
     */
    final private static List<String> tabCompletionTestPackages;

    // Static block to initialize tabCompletionTestPackages:
    static {
        final ClassFinder classFinder = new ClassFinder();
        classFinder.find("org.broadinstitute.hellbender", CommandLineProgram.class);
        tabCompletionTestPackages = classFinder.getClasses().stream().map(Class::getName).collect(Collectors.toList());
    }

    @Test
    public static void documentationSmokeTest() throws IOException, InterruptedException {
        final File tabCompletionTestTarget = createTempDir("tabCompletionTest");

        final String[] argArray = new String[] {

                "-doclet", BashTabCompletionDoclet.class.getName(),

                "-docletpath", "build/libs/",
                "-sourcepath", "src/main/java",
                "-d", tabCompletionTestTarget.getAbsolutePath(), // directory must exist
                "-cp", System.getProperty("java.class.path"),

                "-use-default-templates",

                "-output-file-extension", "sh",
                "-index-file-extension", "sh",
                "-absolute-version", "0.0-001",
                "-build-timestamp", new SimpleDateFormat("dd-mm-yyyy hh:mm:ss").format( new Date() ),

                "-caller-script-name", "gatk-launch",

                "-caller-pre-legal-args", "--help --list --dryRun --javaOptions",
                "-caller-pre-arg-val-types", "null null null String",
                "-caller-pre-mutex-args", "--help;list,dryRun,javaOptions --list;help,dryRun,javaOptions",
                "-caller-pre-alias-args", "--help;-h",
                "-caller-pre-arg-min-occurs", "0 0 0 0",
                "-caller-pre-arg-max-occurs", "1 1 1 1",

                "-caller-post-legal-args", "--sparkRunner --spark-master --cluster --dryRun --javaOptions --conf --driver-memory --driver-cores --executor-memory --executor-cores --num-executors",
                "-caller-post-arg-val-types", "String String String null String file int int int int int",
                "-caller-post-mutex-args", "",
                "-caller-post-alias-args", "",
                "-caller-post-arg-min-occurs", "0 0 0 0 0 0 0 0 0 0",
                "-caller-post-arg-max-occurs", "1 1 1 1 1 1 1 1 1 1",

                "-verbose"
        };

        final List<String> docArgList = new ArrayList<>();
        docArgList.addAll(Arrays.asList(argArray));
        docArgList.addAll(tabCompletionTestPackages);

        // This is  smoke test; we just want to make sure it doesn't blow up
        int success = com.sun.tools.javadoc.Main.execute(docArgList.toArray(new String[]{}));
        Assert.assertEquals(success, 0, "Failure processing gatkTabComplete via javadoc");
    }

    @Test
    public static void documentationSmokeTestFixed() throws IOException, InterruptedException {
        final File tabCompletionTestTarget = createTempDir("tabCompletionTest");

        final List<String> argList = Arrays.asList(

                "javadoc",

                "-doclet", BashTabCompletionDoclet.class.getName(),

                "-docletpath", "build/libs/",
                "-sourcepath", "src/main/java",
                "-d", tabCompletionTestTarget.getAbsolutePath(), // directory must exist
                "-cp", System.getProperty("java.class.path"),

                "-use-default-templates",

                "-output-file-extension", "sh",
                "-index-file-extension", "sh",
                "-absolute-version", "0.0-001",
                "-build-timestamp", new SimpleDateFormat("dd-mm-yyyy hh:mm:ss").format( new Date() ),

                "-caller-script-name", "gatk-launch",

                "-caller-pre-legal-args", "--help --list --dryRun --javaOptions",
                "-caller-pre-arg-val-types", "null null null String",
                "-caller-pre-mutex-args", "--help;list,dryRun,javaOptions --list;help,dryRun,javaOptions",
                "-caller-pre-alias-args", "--help;-h",
                "-caller-pre-arg-min-occurs", "0 0 0 0",
                "-caller-pre-arg-max-occurs", "1 1 1 1",

                "-caller-post-legal-args", "--sparkRunner --spark-master --cluster --dryRun --javaOptions --conf --driver-memory --driver-cores --executor-memory --executor-cores --num-executors",
                "-caller-post-arg-val-types", "String String String null String file int int int int int",
                "-caller-post-mutex-args", "",
                "-caller-post-alias-args", "",
                "-caller-post-arg-min-occurs", "0 0 0 0 0 0 0 0 0 0",
                "-caller-post-arg-max-occurs", "1 1 1 1 1 1 1 1 1 1",

                "-verbose"
        );

        // This is  smoke test; we just want to make sure it doesn't blow up

        // Run this as a process, not through Java itself:
        final ProcessController processController = new ProcessController();
        runProcess(processController, argList.toArray(new String[] {}) );
    }

}
