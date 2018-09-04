package org.broadinstitute.hellbender.tools.examples;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.tools.examples.metrics.multi.ExampleCollectMultiMetricsSpark;
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
public final class ExampleCollectMultiMetricsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "/examples/metrics");

    @Override
    public String getTestedClassName() {
        return ExampleCollectMultiMetricsSpark.class.getSimpleName();
    }

    @DataProvider(name="metricsInputFiles")
    public Object[][] insertSizeMetricsInputFiles() {
        return new Object[][] {
                {"exampleMetrics.bam", true, "expectedMultiMetricsL3.txt"},
                {"exampleMetrics.bam", false, "expectedMultiMetricsL1.txt"},
        };
    }
    @Test(dataProvider="metricsInputFiles", groups = "spark")
    public void test(
            final String fileName,
            final boolean allLevels, // collect metrics for each possible level
            final String expectedResultsFile) throws IOException {

        final File input = new File(TEST_DATA_DIR, fileName);
        final File textOut = GATKBaseTest.createTempFile("test", ".txt");

        final ArgumentsBuilder args = new ArgumentsBuilder();

        // IO arguments
        args.addInput(input);
        args.addOutput(textOut);

        if (allLevels) {
            args.addArgument("LEVEL", MetricAccumulationLevel.ALL_READS.toString());
            args.addArgument("LEVEL", MetricAccumulationLevel.SAMPLE.toString());
            args.addArgument("LEVEL", MetricAccumulationLevel.LIBRARY.toString());
            args.addArgument("LEVEL", MetricAccumulationLevel.READ_GROUP.toString());
        }

        this.runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(
                textOut,
                new File(TEST_DATA_DIR, expectedResultsFile),
                "#"
        );
    }
}
