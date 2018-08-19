package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.examples.metrics.single.ExampleCollectSingleMetricsSpark;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Test for example multi-level metrics collector.
 */
public final class ExampleCollectSingleMetricsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "/examples/metrics");

    @Override
    public String getTestedClassName() {
        return ExampleCollectSingleMetricsSpark.class.getSimpleName();
    }

    @DataProvider(name="metricsInputFiles")
    public Object[][] insertSizeMetricsInputFiles() {
        return new Object[][] {
                {"exampleMetrics.bam", "expectedSingleMetrics.txt"},
        };
    }
    @Test(dataProvider="metricsInputFiles", groups = "spark")
    public void test(
            final String fileName,
            final String expectedResultsFile) throws IOException {

        final File input = new File(TEST_DATA_DIR, fileName);
        final File textOut = GATKBaseTest.createTempFile("test", ".txt");

        final ArgumentsBuilder args = new ArgumentsBuilder();

        // IO arguments
        args.addInput(input);
        args.addOutput(textOut);

        this.runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(
                textOut,
                new File(TEST_DATA_DIR, expectedResultsFile),
                "#"
        );
    }
}
