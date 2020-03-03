package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;

/**
 * Integration tests for {@link CollectInsertSizeMetricsSpark}.
 */
public final class CollectInsertSizeMetricsSparkIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(
            "src/test/resources/org/broadinstitute/hellbender/metrics/analysis/CollectInsertSizeMetrics");

    @Override
    public String getTestedClassName() {
        return CollectInsertSizeMetricsSpark.class.getSimpleName();
    }

    @DataProvider(name="metricsInputFiles")
    public Object[][] insertSizeMetricsInputFiles() {
        return new Object[][] {
                // single level collection
                {"insert_size_metrics_test.sam", null, false, "expectedInsertSizeMetricsL1.txt"},
                {"insert_size_metrics_test.bam", null, false, "expectedInsertSizeMetricsL1.txt"},
                {"insert_size_metrics_test.cram", hg19_chr1_1M_Reference, false, "expectedInsertSizeMetricsL1.txt"},

                // collect for all levels
                {"insert_size_metrics_test.sam", null, true, "expectedInsertSizeMetricsL3.txt"},
                {"insert_size_metrics_test.bam", null, true, "expectedInsertSizeMetricsL3.txt"},
                {"insert_size_metrics_test.cram", hg19_chr1_1M_Reference, true, "expectedInsertSizeMetricsL3.txt"}
        };
    }

    //NOTE: These tests run on small files in a test environment, so only a single partition is
    // used, and the collector's combine/reduce code is never executed. In order to test that,
    // these same tests are included in in InsertSizeMetricsCollectorSparkUnitTest, using the
    // InsertSizeMetricsCollectorSpark object directly, where a repartition on the RDD is used to
    // force the code to run on more than one partition and execute the combine/reduce code.
    @Test(dataProvider="metricsInputFiles", groups = {"R", "spark"})
    public void test(
            final String fileName,
            final String referenceName,
            final boolean allLevels, // collect metrics for each possible level
            final String expectedResultsFile) throws IOException {

        // set up test data input and result outputs (two: one text one histogram plot in pdf)
        final File input = new File(TEST_DATA_DIR, fileName);

        final File textOut = GATKBaseTest.createTempFile("test", ".txt");
        final File pdfOut = GATKBaseTest.createTempFile("test", ".pdf");

        final ArgumentsBuilder args = new ArgumentsBuilder();

        // IO arguments
        args.addRaw("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        args.addRaw(input.getAbsolutePath());

        args.addRaw("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        args.addRaw(textOut.getAbsolutePath());

        args.addRaw("--produce-plot");
        args.addRaw("-" + "histogram-plot-file");
        args.addRaw(pdfOut.getAbsolutePath());

        if (null != referenceName) {
            final File REF = new File(referenceName);
            args.addRaw("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
            args.addRaw(REF.getAbsolutePath());
        }

        // some filter options
        args.addRaw("-DF");
        args.addRaw(ReadFilterLibrary.FirstOfPairReadFilter.class.getSimpleName());
        args.addRaw("-RF");
        args.addRaw(ReadFilterLibrary.SecondOfPairReadFilter.class.getSimpleName());

        if (allLevels) {
            // accumulation level options (all included for better test coverage)
            args.addRaw("-" + "LEVEL");
            args.addRaw(MetricAccumulationLevel.ALL_READS.toString());
            args.addRaw("-" + "LEVEL");
            args.addRaw(MetricAccumulationLevel.SAMPLE.toString());
            args.addRaw("-" + "LEVEL");
            args.addRaw(MetricAccumulationLevel.LIBRARY.toString());
            args.addRaw("-" + "LEVEL");
            args.addRaw(MetricAccumulationLevel.READ_GROUP.toString());
        }

        this.runCommandLine(args.getArgsArray());

        IntegrationTestSpec.assertEqualTextFiles(
                textOut,
                new File(TEST_DATA_DIR, expectedResultsFile),
                "#"
        );
    }
}
