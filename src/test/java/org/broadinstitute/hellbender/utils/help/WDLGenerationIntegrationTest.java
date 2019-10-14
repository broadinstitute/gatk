package org.broadinstitute.hellbender.utils.help;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Smoke test to run doc gen on a subset of classes to make sure it doesn't regress.
 */
public class WDLGenerationIntegrationTest extends CommandLineProgramTest {
    /**
     * Entry point for manually running the gatkWDLGen process on a subset of packages from within GATK.
     */
    private static String[] docTestPackages = {
            "org.broadinstitute.hellbender.cmdline.argumentcollections",
            "org.broadinstitute.hellbender.cmdline.GATKPlugin",
            "org.broadinstitute.hellbender.engine.filters",
            "org.broadinstitute.hellbender.tools",
            "org.broadinstitute.hellbender.tools.spark.transforms.bqsr",
            "org.broadinstitute.hellbender.tools.spark.transforms.markduplicates",
            "org.broadinstitute.hellbender.tools.walkers.bqsr",
            "org.broadinstitute.hellbender.tools.walkers.vqsr",
            "org.broadinstitute.hellbender.tools.walkers.variantutils",
            //not including any Picard tools yet
    };

    @Test
    public static void wdlGenSmokeTest() {
        File docTestTarget = createTempDir("wdlgentest");
        String[] argArray = new String[]{
                "javadoc",
                "-doclet", GATKWDLDoclet.class.getName(),
                "-docletpath", System.getProperty("java.class.path"),
                "-sourcepath", "src/main/java",
                "-settings-dir", "src/main/resources/org/broadinstitute/hellbender/utils/wdlTemplates",
                "-d", docTestTarget.getAbsolutePath(), // directory must exist
                "-output-file-extension", "wdl",
                "-build-timestamp", "2016/11/11 11:11:11",
                "-absolute-version", "1.1-111",
                "-cp", System.getProperty("java.class.path"),
                "-verbose"
        };

        final List<String> docArgList = new ArrayList<>();
        docArgList.addAll(Arrays.asList(argArray));
        docArgList.addAll(Arrays.asList(docTestPackages));

        // This is  smoke test; we just want to make sure it doesn't blow up

        // Run this as a process, not through Java itself:
        final ProcessController processController = new ProcessController();
        runProcess(processController, docArgList.toArray(new String[] {}), "Failure processing gatkWDL via javadoc" );
    }
}
