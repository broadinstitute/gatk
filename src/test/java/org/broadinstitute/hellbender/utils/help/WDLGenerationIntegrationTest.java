package org.broadinstitute.hellbender.utils.help;

import org.broadinstitute.barclay.argparser.ClassFinder;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;

/**
 * Smoke test to run doc gen on a subset of classes to make sure it doesn't regress.
 */
public class WDLGenerationIntegrationTest extends CommandLineProgramTest {

    final private static List<String> wdlGenTestPackages;
    static {
        final ClassFinder classFinder = new ClassFinder();
        classFinder.find("org.broadinstitute.hellbender", CommandLineProgram.class);
        wdlGenTestPackages = Collections.unmodifiableList(
                classFinder.getClasses()
                        .stream()
                        .map(cl -> cl.getPackage().getName())
                        .collect(Collectors.toSet()) // uniquify
                        .stream()
                        .collect(Collectors.toList())
        );
    }

    @Test
    public static void wdlGenSmokeTest() {
        final File wdlTestTargetDir = createTempDir("wdlgentest");
        doWDLGenTest(wdlGenTestPackages, "src/main/java", wdlTestTargetDir);
    }

    @Test
    public static void wdlGenTemplateTest() throws IOException {
        final File expectedResultsDir = new File("src/test/resources/org/broadinstitute/hellbender/utils/wdltest/");
        final File wdlTestTargetDir = createTempDir("wdlgentemplatetest");

        doWDLGenTest(
                Collections.singletonList("org.broadinstitute.hellbender.utils.help"),
                "src/test/java",
                wdlTestTargetDir);

        IntegrationTestSpec.assertEqualTextFiles(
                new File(expectedResultsDir, "index.html"),
                new File(wdlTestTargetDir, "index.html"));
        IntegrationTestSpec.assertEqualTextFiles(
                new File(expectedResultsDir, "TestWDLTool.wdl"),
                new File(wdlTestTargetDir, "TestWDLTool.wdl"));
        IntegrationTestSpec.assertEqualTextFiles(
                new File(expectedResultsDir, "TestWDLToolAllArgs.wdl"),
                new File(wdlTestTargetDir, "TestWDLToolAllArgs.wdl"));
        IntegrationTestSpec.assertEqualTextFiles(
                new File(expectedResultsDir, "TestWDLToolAllArgsTest.wdl"),
                new File(wdlTestTargetDir, "TestWDLToolAllArgsTest.wdl"));
        IntegrationTestSpec.assertEqualTextFiles(
                new File(expectedResultsDir, "TestWDLToolInputs.json"),
                new File(wdlTestTargetDir, "TestWDLToolInputs.json"));
        IntegrationTestSpec.assertEqualTextFiles(
                new File(expectedResultsDir, "TestWDLToolAllArgsInputs.json"),
                new File(wdlTestTargetDir, "TestWDLToolAllArgsInputs.json"));
        IntegrationTestSpec.assertEqualTextFiles(
                new File(expectedResultsDir, "TestWDLToolAllArgsTestInputs.json"),
                new File(wdlTestTargetDir, "TestWDLToolAllArgsTestInputs.json"));
    }

    // suppress deprecation warning on Java 11 since we're using deprecated javadoc APIs
    @SuppressWarnings({"deprecation","removal"})
    public static void doWDLGenTest(List<String> testPackages, final String sourcePath, final File wdlTestTargetDir) {

        final String[] argArray = new String[]{
                "javadoc",
                "-doclet", GATKWDLDoclet.class.getName(),
                "-docletpath", System.getProperty("java.class.path"),
                "-sourcepath", sourcePath,
                "-settings-dir", "src/main/resources/org/broadinstitute/hellbender/utils/wdlTemplates",
                "-d", wdlTestTargetDir.getAbsolutePath(), // directory must exist
                "-output-file-extension", "wdl",
                "-build-timestamp", "2016/11/11 11:11:11",
                "-build-dir", ".",
                "-absolute-version", "1.1-111",
                "-verbose"
        };

        final List<String> docArgList = new ArrayList<>();
        docArgList.addAll(Arrays.asList(argArray));
        docArgList.add("-cp");
        docArgList.add(System.getProperty("java.class.path"));
        docArgList.addAll(testPackages);

        // Run javadoc in the current JVM with the custom WDL doclet. This is a smoke test; we just want to
        // make sure it doesn't blow up (the gradle task gatkWDLGenValidation does womtool validation on the results).
        Assert.assertEquals(com.sun.tools.javadoc.Main.execute(docArgList.toArray(new String[] {})), 0);
    }

}
