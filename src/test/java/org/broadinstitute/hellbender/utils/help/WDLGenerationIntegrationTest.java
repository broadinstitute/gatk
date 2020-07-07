package org.broadinstitute.hellbender.utils.help;

import org.broadinstitute.barclay.argparser.ClassFinder;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

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

    // suppress deprecation warning on Java 11 since we're using deprecated javadoc APIs
    @SuppressWarnings({"deprecation","removal"})
    @Test
    public static void wdlGenSmokeTest() {
        final File wdlTestTarget = createTempDir("wdlgentest");

        final String[] argArray = new String[]{
                "javadoc",
                "-doclet", GATKWDLDoclet.class.getName(),
                "-docletpath", System.getProperty("java.class.path"),
                "-sourcepath", "src/main/java",
                "-settings-dir", "src/main/resources/org/broadinstitute/hellbender/utils/wdlTemplates",
                "-d", wdlTestTarget.getAbsolutePath(), // directory must exist
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
        docArgList.addAll(wdlGenTestPackages);

        // Run javadoc in the current JVM with the custom WDL doclet. This is a smoke test; we just want to
        // make sure it doesn't blow up (the gradle task gatkWDLGenValidation does womtool validation on the results).
        Assert.assertEquals(com.sun.tools.javadoc.Main.execute(docArgList.toArray(new String[] {})), 0);
    }
}
