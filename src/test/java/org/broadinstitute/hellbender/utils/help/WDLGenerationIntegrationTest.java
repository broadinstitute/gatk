//package org.broadinstitute.hellbender.utils.help;
//
//import org.json.simple.parser.JSONParser;
//import org.broadinstitute.barclay.argparser.ClassFinder;
//import org.broadinstitute.hellbender.CommandLineProgramTest;
//import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
//import org.json.simple.parser.ParseException;
//import org.testng.Assert;
//import org.testng.annotations.Test;
//
//import java.io.File;
//import java.io.FileReader;
//import java.io.IOException;
//import java.util.ArrayList;
//import java.util.Arrays;
//import java.util.Collections;
//import java.util.List;
//import java.util.stream.Collectors;
//
//import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
//
///**
// * Smoke test to run doc gen on a subset of classes to make sure it doesn't regress.
// */
//public class WDLGenerationIntegrationTest extends CommandLineProgramTest {
//
//    final private static List<String> wdlGenTestPackages;
//    static {
//        final ClassFinder classFinder = new ClassFinder();
//        classFinder.find("org.broadinstitute.hellbender", CommandLineProgram.class);
//        wdlGenTestPackages = Collections.unmodifiableList(
//                classFinder.getClasses()
//                        .stream()
//                        .map(cl -> cl.getPackage().getName())
//                        .collect(Collectors.toSet()) // uniquify
//                        .stream()
//                        .collect(Collectors.toList())
//        );
//    }
//
//    @Test
//    private void wdlGenSmokeTest() throws IOException, ParseException {
//        final File wdlTestTargetDir = createTempDir("wdlgentest");
//        doWDLGenTest(wdlGenTestPackages, "src/main/java", wdlTestTargetDir);
//
//        // load and parse every generated JSON file to make sure they're valid JSON
//        final File[] jsonFiles = wdlTestTargetDir.listFiles((File dir, String name) -> name.endsWith(".json"));
//        for (final File f : jsonFiles) {
//            assertValidJSONFile(f);
//        }
//    }
//
//    // This test uses a test tool with all combinations of input/output, list/scalar, across all types, including
//    // companions, to ensure that the various templates use the template maps correctly.
//    @Test
//    private void wdlGenTemplateTest() throws IOException, ParseException {
//        final File expectedResultsDir = new File("src/test/resources/org/broadinstitute/hellbender/utils/WDLTestExpectedOutputs/");
//        final File wdlTestTargetDir = createTempDir("wdlgentemplatetest");
//
//        doWDLGenTest(
//                Collections.singletonList("org.broadinstitute.hellbender.utils.help"),
//                "src/test/java",
//                wdlTestTargetDir);
//
//        // index
//        final String indexHTML = "index.html";
//        IntegrationTestSpec.assertEqualTextFiles(
//                new File(expectedResultsDir, "expected" + indexHTML),
//                new File(wdlTestTargetDir, indexHTML));
//
//        // wdls
//        final String defaultWDL = "TestWDLTool.wdl";
//        IntegrationTestSpec.assertEqualTextFiles(
//                new File(expectedResultsDir, "expected" + defaultWDL),
//                new File(wdlTestTargetDir, defaultWDL));
//
//        final String allArgsWDL = "TestWDLToolAllArgs.wdl";
//        IntegrationTestSpec.assertEqualTextFiles(
//                new File(expectedResultsDir, "expected" + allArgsWDL),
//                new File(wdlTestTargetDir, allArgsWDL));
//
//        final String allArgsTestWDL = "TestWDLToolAllArgsTest.wdl";
//        IntegrationTestSpec.assertEqualTextFiles(
//                new File(expectedResultsDir, "expected" + allArgsTestWDL),
//                new File(wdlTestTargetDir, allArgsTestWDL));
//
//        // jsons
//        final String defaultWDLInputs = "TestWDLToolInputs.json";
//        IntegrationTestSpec.assertEqualTextFiles(
//                new File(expectedResultsDir, "expected" + defaultWDLInputs),
//                new File(wdlTestTargetDir, defaultWDLInputs));
//        assertValidJSONFile(new File(wdlTestTargetDir, defaultWDLInputs));
//
//        final String allArgsWDLInputs = "TestWDLToolAllArgsInputs.json";
//        IntegrationTestSpec.assertEqualTextFiles(
//                new File(expectedResultsDir, "expected" + allArgsWDLInputs),
//                new File(wdlTestTargetDir, allArgsWDLInputs));
//        assertValidJSONFile(new File(wdlTestTargetDir, allArgsWDLInputs));
//
//        final String allArgsTestWDLInputs = "TestWDLToolAllArgsTestInputs.json";
//        IntegrationTestSpec.assertEqualTextFiles(
//                new File(expectedResultsDir, "expected" + allArgsTestWDLInputs),
//                new File(wdlTestTargetDir, allArgsTestWDLInputs));
//        assertValidJSONFile(new File(wdlTestTargetDir, allArgsTestWDLInputs));
//    }
//
//    private void assertValidJSONFile(final File targetFile) throws IOException, ParseException {
//        try (FileReader fileReader = new FileReader(targetFile)) {
//            new JSONParser().parse(fileReader);
//        }
//    }
//
//    // suppress deprecation warning on Java 11 since we're using deprecated javadoc APIs
//    @SuppressWarnings({"deprecation","removal"})
//    private void doWDLGenTest(List<String> testPackages, final String sourcePath, final File wdlTestTargetDir) {
//
//        final String[] argArray = new String[]{
//                "javadoc",
//                "-doclet", GATKWDLDoclet.class.getName(),
//                "-docletpath", System.getProperty("java.class.path"),
//                "-sourcepath", sourcePath,
//                "-settings-dir", "src/main/resources/org/broadinstitute/hellbender/utils/wdlTemplates",
//                "-d", wdlTestTargetDir.getAbsolutePath(), // directory must exist
//                "-output-file-extension", "wdl",
//                "-build-timestamp", "2016/11/11 11:11:11",
//                "-build-dir", ".",
//                "-absolute-version", "1.2.3.4",
//                "-verbose"
//        };
//
//        final List<String> docArgList = new ArrayList<>();
//        docArgList.addAll(Arrays.asList(argArray));
//        docArgList.add("-cp");
//        docArgList.add(System.getProperty("java.class.path"));
//        docArgList.addAll(testPackages);
//
//        // Run javadoc in the current JVM with the custom WDL doclet. This is a smoke test; we just want to
//        // make sure it doesn't blow up (the gradle task gatkWDLGenValidation does womtool validation on the results).
//        Assert.assertEquals(com.sun.tools.javadoc.Main.execute(docArgList.toArray(new String[] {})), 0);
//    }
//
//}
