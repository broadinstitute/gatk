package org.broadinstitute.hellbender.utils.help;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.Test;

import java.util.spi.ToolProvider;
import java.io.File;
import java.io.IOException;
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
            "org.broadinstitute.hellbender.cmdline.argumentcollections",
            "org.broadinstitute.hellbender.cmdline.GATKPlugin",
            "org.broadinstitute.hellbender.engine.filters",
            "org.broadinstitute.hellbender.tools",
            "org.broadinstitute.hellbender.tools.spark",
            "org.broadinstitute.hellbender.tools.spark.pipelines",
            "org.broadinstitute.hellbender.tools.spark.pipelines.metrics",
            //"org.broadinstitute.hellbender.tools.spark.transforms.bqsr",
            "org.broadinstitute.hellbender.tools.spark.transforms.markduplicates",
            "org.broadinstitute.hellbender.tools.walkers.bqsr",
            "org.broadinstitute.hellbender.tools.walkers.vqsr",
            "org.broadinstitute.hellbender.tools.walkers.variantutils",
            //"picard.fingerprint",
            //"picard.analysis"
    };

    @Test
    public void documentationSmokeTest() throws IOException {
        if (!isGATKDockerContainer()) {
            final File docTestTarget = createTempDir("docgentest");
            final String[] argArray = new String[]{
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
            docArgList.addAll(Arrays.asList(docTestPackages));

            // Run javadoc in the current JVM with the custom doclet, and make sure it succeeds (this is a smoke test;
            // we just want to make sure it doesn't blow up).
            final ToolProvider jdProvider = ToolProvider.findFirst("javadoc")
                    .orElseThrow(() -> new IllegalStateException("Can't find javadoc tool"));

            try (final StringWriter stringWriter = new StringWriter();
                 final PrintWriter outputWriter = new PrintWriter(stringWriter)) {

                final String[] args = docArgList.toArray(new String[]{});
                final int retCode = jdProvider.run(outputWriter, outputWriter, args);

                // just make sure the task succeeded
                Assert.assertEquals(retCode, 0, stringWriter.toString());
            }
        }
    }
}
