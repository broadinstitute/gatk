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
public class DocumentationGenerationIntegrationTest extends CommandLineProgramTest {
    /**
     * Entry point for manually running the gatkDoc process on a subset of packages from within GATK.
     */
    private static String[] docTestPackages = {
            "org.broadinstitute.hellbender.cmdline.argumentcollections",
            "org.broadinstitute.hellbender.cmdline.GATKPlugin",
            "org.broadinstitute.hellbender.engine.filters",
            "org.broadinstitute.hellbender.tools",
            "org.broadinstitute.hellbender.tools.spark",
            "org.broadinstitute.hellbender.tools.spark.pipelines",
            "org.broadinstitute.hellbender.tools.spark.pipelines.metrics",
            "org.broadinstitute.hellbender.tools.spark.transforms.bqsr",
            "org.broadinstitute.hellbender.tools.spark.transforms.markduplicates",
            "org.broadinstitute.hellbender.tools.walkers.bqsr",
            "org.broadinstitute.hellbender.tools.walkers.vqsr",
            "org.broadinstitute.hellbender.tools.walkers.variantutils",
            "picard.fingerprint",
            "picard.analysis"
    };

    @Test
    public static void documentationSmokeTest() throws IOException {
        File docTestTarget = createTempDir("docgentest");
        String[] argArray = new String[]{
                "javadoc",
                "-doclet", GATKHelpDoclet.class.getName(),
                "-docletpath", System.getProperty("java.class.path"),
                "-sourcepath", "src/main/java",
                "-settings-dir", "src/main/resources/org/broadinstitute/hellbender/utils/helpTemplates",
                "-d", docTestTarget.getAbsolutePath(), // directory must exist
                "-output-file-extension", "html",
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
        runProcess(processController, docArgList.toArray(new String[] {}), "Failure processing gatkDoc via javadoc" );
    }
}
