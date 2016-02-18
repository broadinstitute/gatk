package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetrics;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public final class CollectInsertSizeMetricsSparkUnitTest extends CommandLineProgramTest{

    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/analysis/CollectInsertSizeMetrics");

    public String getTestedClassName() {
        return CollectInsertSizeMetricsSpark.class.getSimpleName();
    }

    @DataProvider(name="metricsfiles")
    public Object[][] insertSizeMetricsFiles() {
        return new Object[][] {
                {"insert_size_metrics_test.bam", null},
                {"insert_size_metrics_test.cram", hg19_chr1_1M_Reference} // TODO: cram requires ref, but actually ref not quite used yet currently
        };
    }

    @Test(dataProvider="metricsfiles")
    public void test(final String fileName, final String referenceName) throws IOException {

        // set up test data input and result outputs (two: one text one histogram plot in pdf)
        final File input = new File(TEST_DATA_DIR, fileName);
        final File textOut = BaseTest.createTempFile("test", ".txt");
        final File pdfOut = BaseTest.createTempFile("test", ".pdf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        // IO arguments
        args.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        args.add(input.getAbsolutePath());
        args.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        args.add(textOut.getAbsolutePath());
        args.add("-" + "HIST");
        args.add(pdfOut.getAbsolutePath());
        if (null != referenceName) {
            final File REF = new File(referenceName);
            args.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
            args.add(REF.getAbsolutePath());
        }

        // some filter options
        args.add("-" + "E");
        args.add(CollectInsertSizeMetricsSpark.UseEnd.SECOND);

        this.runCommandLine(args.getArgsArray());

        final MetricsFile<InsertSizeMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(textOut));

        Assert.assertFalse(output.getMetrics().isEmpty());

        final double DOUBLE_TOLERANCE = 0.05;

        for (final InsertSizeMetrics metrics : output.getMetrics()) {

            // TODO: not yet checked
            // Assert.assertEquals(metrics.PAIR_ORIENTATION.name(), "FR");

            // TODO: add tests for different collection level
            if (metrics.LIBRARY == null) {  // SAMPLE or ALL_READS level

                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 45);
                Assert.assertEquals(metrics.READ_PAIRS, 13);

                Assert.assertEquals(metrics.MEAN_INSERT_SIZE, 40.1, DOUBLE_TOLERANCE);
                Assert.assertEquals(metrics.STANDARD_DEVIATION, 3.12, DOUBLE_TOLERANCE);

                Assert.assertEquals((int) metrics.MEDIAN_INSERT_SIZE, 41);
                Assert.assertEquals((int) metrics.MEDIAN_ABSOLUTE_DEVIATION, 3);

                // following is also a test on symmetric bin width collection method
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 9);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 11);

            } else {
                Assert.fail("Unexpected metric: " + metrics);
            }
        }
    }
}
