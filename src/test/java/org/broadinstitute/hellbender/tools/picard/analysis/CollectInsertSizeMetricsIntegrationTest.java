package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.metrics.InsertSizeMetrics;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * Tests multi-level CollectInsertSizeMetrics
 */
public final class CollectInsertSizeMetricsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/analysis/CollectInsertSizeMetrics");

    @Override
    public String getTestedClassName() {
        return CollectInsertSizeMetrics.class.getSimpleName();
    }

    @DataProvider(name="metricsfiles")
    public Object[][] insertSizeMetricsFiles() {
        return new Object[][] {
                {"insert_size_metrics_test.sam", null, false, "expectedInsertSizeMetricsL1.txt"},
                {"insert_size_metrics_test.bam", null, false, "expectedInsertSizeMetricsL1.txt"},
                {"insert_size_metrics_test.cram", hg19_chr1_1M_Reference, false, "expectedInsertSizeMetricsL1.txt"},

                {"insert_size_metrics_test.sam", null, true, "expectedInsertSizeMetricsL3.txt"},
                {"insert_size_metrics_test.bam", null, true, "expectedInsertSizeMetricsL3.txt"},
                {"insert_size_metrics_test.cram", hg19_chr1_1M_Reference, true, "expectedInsertSizeMetricsL3.txt"}
        };
    }

    @Test(dataProvider="metricsfiles")
    public void test(
            final String fileName,
            final String referenceName,
            final boolean allLevels,
            final String expectedResultsFile) throws IOException {
        final File input = new File(TEST_DATA_DIR, fileName);
        final File outfile = BaseTest.createTempFile("test", ".insert_size_metrics");
        final File pdf = BaseTest.createTempFile("test", ".pdf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(input.getAbsolutePath());
        args.add("--output");
        args.add(outfile.getAbsolutePath());
        if (null != referenceName) {
            final File REF = new File(referenceName);
            args.add("-R");
            args.add(REF.getAbsolutePath());
        }

        if (allLevels) {
            // accumulation level options (all included for better test coverage)
            args.add("-" + "LEVEL");
            args.add(MetricAccumulationLevel.ALL_READS.toString());
            args.add("-" + "LEVEL");
            args.add(MetricAccumulationLevel.SAMPLE.toString());
            args.add("-" + "LEVEL");
            args.add(MetricAccumulationLevel.LIBRARY.toString());
            args.add("-" + "LEVEL");
            args.add(MetricAccumulationLevel.READ_GROUP.toString());
        }

        args.add("--producePlot");
        args.add("--histogramPlotFile");
        args.add(pdf.getAbsolutePath());

        runCommandLine(args.getArgsArray());

        IntegrationTestSpec.assertEqualTextFiles(
                outfile,
                new File(TEST_DATA_DIR, expectedResultsFile),
                "#"
        );
    }
}
