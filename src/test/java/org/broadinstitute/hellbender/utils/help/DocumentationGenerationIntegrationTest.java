package org.broadinstitute.hellbender.utils.help;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Smoke test to run doc gen on a subset of classes to make sure it doesn't regress.
 */
public class DocumentationGenerationIntegrationTest extends CommandLineProgramTest {
    /**
     * Entry point for manually running the gatkDoc process on a subset of packages from within GATK.
     */
    private static String[] docTestPackages = {
            //"org.broadinstitute.hellbender.cmdline.argumentcollections",
            //"org.broadinstitute.hellbender.cmdline.GATKPlugin",
            //"org.broadinstitute.hellbender.engine.filters",
            //"org.broadinstitute.hellbender.tools",
            //"org.broadinstitute.hellbender.tools.spark",
            //"org.broadinstitute.hellbender.tools.spark.pipelines",
            //"org.broadinstitute.hellbender.tools.spark.pipelines.metrics",
            // binary search: the error is somewhere below here
            //"org.broadinstitute.hellbender.tools.spark.transforms.bqsr",
            //"org.broadinstitute.hellbender.tools.spark.transforms.markduplicates",
            //"org.broadinstitute.hellbender.tools.walkers.bqsr",
            // more binary search: the error is somewhere below here
            "org.broadinstitute.hellbender.tools.walkers.vqsr",
            "org.broadinstitute.hellbender.tools.walkers.variantutils",
            "picard.fingerprint",
            "picard.analysis"
    };

    private static String[] docTestPackages1 = {
            "org.broadinstitute.hellbender.tools.walkers.vqsr"
            // "org.broadinstitute.hellbender.tools.walkers.variantutils",
            // "picard.fingerprint",
            // "picard.analysis"
    };

    private static String[] docTestPackages2 = {
            //"org.broadinstitute.hellbender.tools.walkers.vqsr"
            "org.broadinstitute.hellbender.tools.walkers.variantutils"
            // "picard.fingerprint",
            // "picard.analysis"
    };

    private static String[] docTestPackages3 = {
            // "org.broadinstitute.hellbender.tools.walkers.vqsr"
            // "org.broadinstitute.hellbender.tools.walkers.variantutils",
            "picard.fingerprint"
            // "picard.analysis"
    };

    private static String[] docTestPackages4 = {
            // "org.broadinstitute.hellbender.tools.walkers.vqsr"
            // "org.broadinstitute.hellbender.tools.walkers.variantutils",
            // "picard.fingerprint",
            "picard.analysis"
    };

    // suppress deprecation warning on Java 11 since we're using deprecated javadoc APIs
    @SuppressWarnings({"deprecation","removal"})
    @Test
    public static void documentationSmokeTest1() {
        final File docTestTarget = createTempDir("docgentest");
        final String[] argArray = new String[]{
                "javadoc",
                "-doclet", GATKHelpDoclet.class.getName(),
                "-docletpath", System.getProperty("java.class.path"),
                "-sourcepath", "src/main/java",
                "-settings-dir", "src/main/resources/org/broadinstitute/hellbender/utils/helpTemplates",
                "-d", docTestTarget.getAbsolutePath(), // directory must exist
                "-output-file-extension", "html",
                "-build-timestamp", "2016/11/11 11:11:11",
                "-absolute-version", "1.1-111",
                "-verbose"
        };

        final List<String> docArgList = new ArrayList<>();
        docArgList.addAll(Arrays.asList(argArray));
        docArgList.add("-cp");
        docArgList.add(System.getProperty("java.class.path"));
        docArgList.addAll(Arrays.asList(docTestPackages1));

        // Run javadoc in the current JVM with the custom doclet, and make sure it succeeds (this is a smoke test;
        // we just want to make sure it doesn't blow up).
        Assert.assertEquals(com.sun.tools.javadoc.Main.execute(docArgList.toArray(new String[] {})), 0);
    }

    // suppress deprecation warning on Java 11 since we're using deprecated javadoc APIs
    @SuppressWarnings({"deprecation","removal"})
    @Test
    public static void documentationSmokeTest2() {
        final File docTestTarget = createTempDir("docgentest");
        final String[] argArray = new String[]{
                "javadoc",
                "-doclet", GATKHelpDoclet.class.getName(),
                "-docletpath", System.getProperty("java.class.path"),
                "-sourcepath", "src/main/java",
                "-settings-dir", "src/main/resources/org/broadinstitute/hellbender/utils/helpTemplates",
                "-d", docTestTarget.getAbsolutePath(), // directory must exist
                "-output-file-extension", "html",
                "-build-timestamp", "2016/11/11 11:11:11",
                "-absolute-version", "1.1-111",
                "-verbose"
        };

        final List<String> docArgList = new ArrayList<>();
        docArgList.addAll(Arrays.asList(argArray));
        docArgList.add("-cp");
        docArgList.add(System.getProperty("java.class.path"));
        docArgList.addAll(Arrays.asList(docTestPackages2));

        final StringWriter out    = new StringWriter();
        final PrintWriter err = new PrintWriter(out);

        // Run javadoc in the current JVM with the custom doclet, and make sure it succeeds (this is a smoke test;
        // we just want to make sure it doesn't blow up).
        // Assert.assertEquals(com.sun.tools.javadoc.Main.execute(docArgList.toArray(new String[] {})), 0);
        final int result = com.sun.tools.javadoc.Main.execute("programName",
                err, err, err, "docletName",docArgList.toArray(new String[] {}));
        err.flush(); // flush is really optional here, as Writer calls the empty StringWriter.flush
        String message = out.toString();
        Assert.assertEquals(result, 99, message);
    }

    // suppress deprecation warning on Java 11 since we're using deprecated javadoc APIs
    @SuppressWarnings({"deprecation","removal"})
    @Test
    public static void documentationSmokeTest3() {
        final File docTestTarget = createTempDir("docgentest");
        final String[] argArray = new String[]{
                "javadoc",
                "-doclet", GATKHelpDoclet.class.getName(),
                "-docletpath", System.getProperty("java.class.path"),
                "-sourcepath", "src/main/java",
                "-settings-dir", "src/main/resources/org/broadinstitute/hellbender/utils/helpTemplates",
                "-d", docTestTarget.getAbsolutePath(), // directory must exist
                "-output-file-extension", "html",
                "-build-timestamp", "2016/11/11 11:11:11",
                "-absolute-version", "1.1-111",
                "-verbose"
        };

        final List<String> docArgList = new ArrayList<>();
        docArgList.addAll(Arrays.asList(argArray));
        docArgList.add("-cp");
        docArgList.add(System.getProperty("java.class.path"));
        docArgList.addAll(Arrays.asList(docTestPackages3));

        // Run javadoc in the current JVM with the custom doclet, and make sure it succeeds (this is a smoke test;
        // we just want to make sure it doesn't blow up).
        Assert.assertEquals(com.sun.tools.javadoc.Main.execute(docArgList.toArray(new String[] {})), 0);
    }

    // suppress deprecation warning on Java 11 since we're using deprecated javadoc APIs
    @SuppressWarnings({"deprecation","removal"})
    @Test
    public static void documentationSmokeTest4() {
        final File docTestTarget = createTempDir("docgentest");
        final String[] argArray = new String[]{
                "javadoc",
                "-doclet", GATKHelpDoclet.class.getName(),
                "-docletpath", System.getProperty("java.class.path"),
                "-sourcepath", "src/main/java",
                "-settings-dir", "src/main/resources/org/broadinstitute/hellbender/utils/helpTemplates",
                "-d", docTestTarget.getAbsolutePath(), // directory must exist
                "-output-file-extension", "html",
                "-build-timestamp", "2016/11/11 11:11:11",
                "-absolute-version", "1.1-111",
                "-verbose"
        };

        final List<String> docArgList = new ArrayList<>();
        docArgList.addAll(Arrays.asList(argArray));
        docArgList.add("-cp");
        docArgList.add(System.getProperty("java.class.path"));
        docArgList.addAll(Arrays.asList(docTestPackages4));

        // Run javadoc in the current JVM with the custom doclet, and make sure it succeeds (this is a smoke test;
        // we just want to make sure it doesn't blow up).
        Assert.assertEquals(com.sun.tools.javadoc.Main.execute(docArgList.toArray(new String[] {})), 0);
    }
}
