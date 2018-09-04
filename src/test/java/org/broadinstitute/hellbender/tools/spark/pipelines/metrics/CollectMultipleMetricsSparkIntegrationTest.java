package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.Header;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.metrics.MetricsArgumentCollection;
import org.broadinstitute.hellbender.metrics.QualityYieldMetrics;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.metrics.InsertSizeMetrics;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * Integration tests for {@link CollectMultipleMetricsSpark}.
 */
public final class CollectMultipleMetricsSparkIntegrationTest extends CommandLineProgramTest{
    private static final File TEST_DATA_DIR = new File(
            "src/test/resources/org/broadinstitute/hellbender/metrics/analysis/CollectInsertSizeMetrics");

    @Override
    public String getTestedClassName() {
        return CollectMultipleMetricsSpark.class.getSimpleName();
    }

    @DataProvider(name="metricsTestFiles")
    public Object[][] insertSizeMetricsFiles() {
        return new Object[][] {
                // single level collection
                {"insert_size_metrics_test.sam", null, "expectedInsertSizeMetricsL1.txt", "expectedQualityYieldOnInsertSizeMetrics.txt"},
                {"insert_size_metrics_test.bam", null, "expectedInsertSizeMetricsL1.txt", "expectedQualityYieldOnInsertSizeMetrics.txt"},
                {"insert_size_metrics_test.cram", hg19_chr1_1M_Reference, "expectedInsertSizeMetricsL1.txt", "expectedQualityYieldOnInsertSizeMetrics.txt"},
        };
    }

    @Test(dataProvider="metricsTestFiles", groups = "spark")
    public void testBuiltInCollectors(
            final String fileName,
            final String referenceName,
            final String expectedInsertSizeResults,
            final String expectedQualityYieldResults) throws IOException
    {
        ArgumentsBuilder args = new ArgumentsBuilder();
        String outBase = setupMultipleCollector(args, fileName, referenceName);

        // for now, run the only two conforming collectors that we have
        args.add("--collectors" );
        args.add("CollectInsertSizeMetrics" );

        args.add("--collectors" );
        args.add("CollectQualityYieldMetrics" );

        this.runCommandLine(args.getArgsArray());

        validateInsertSizeMetrics(outBase, expectedInsertSizeResults);
        validateQualityYieldMetrics(outBase, expectedQualityYieldResults);
    }

    private String setupMultipleCollector(
            final ArgumentsBuilder args,
            final String fileName,
            final String referenceName) throws IOException
    {
        // set up test data input and result output
        final File input = new File(TEST_DATA_DIR, fileName);

        // create a directory to contain the results since there will be multiple collectors
        // and each may create multiple files
        final File outDir = GATKBaseTest.createTempDir("collectMultiMetricsSparkTest" );
        String outBase = outDir.getAbsolutePath() + "/collectMultiSparkMetrics";

        // IO arguments
        args.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        args.add(input.getAbsolutePath());

        // The output arg specifies only the basename from which output file(s)
        // will derived
        args.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        args.add(outBase);

        if (null != referenceName) {
            final File REF = new File(referenceName);
            args.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
            args.add(REF.getAbsolutePath());
        }

        return outBase;
    }

    private void validateQualityYieldMetrics(final String outBase, final String expectedResults) throws IOException {
        String localOut = outBase + "." + QualityYieldMetrics.getUniqueNameSuffix() + ".txt";

        IntegrationTestSpec.assertEqualTextFiles(
                new File(localOut),
                new File(TEST_DATA_DIR, expectedResults),
                "#");
    }

    private void validateInsertSizeMetrics(final String outBase, final String expectedResults) throws IOException {
        String localOut = outBase + "." + InsertSizeMetrics.getUniqueNameSuffix() + ".txt";

        IntegrationTestSpec.assertEqualTextFiles(
                new File(localOut),
                new File(TEST_DATA_DIR, expectedResults),
                "#");
    }

    // Test implementation of MetricsCollectorSpark used for testing CollectMultipleMetricsSpark
    // with a custom collector added programmatically
    public static class TestCustomCollector implements MetricsCollectorSpark<MetricsArgumentCollection> {
        private static final long serialVersionUID = 1L;
        long count = 0;
        @Override
        public void initialize(
                MetricsArgumentCollection inputArgs, SAMFileHeader samHeader, List<Header> defaultHeaders) {}
        @Override
        public void collectMetrics(JavaRDD<GATKRead> filteredReads, SAMFileHeader samHeader) {
            count = filteredReads.count();
        }
        @Override
        public void saveMetrics(String inputBaseName) {
            //no-op
        }
    }

    @Test(dataProvider="metricsTestFiles", groups = "spark")
    public void testCustomCollectorAPI(
        final String fileName,
        final String referenceName,
        final String expectedInsertSizeResults,
        final String expectedQualityYieldResults) throws IOException
    {
        // Test CollectMultipleMetricsSpark with a custom collector
        final TestCustomCollector testCollector = new TestCustomCollector();

        ArgumentsBuilder args = new ArgumentsBuilder();
        setupMultipleCollector(args, fileName, referenceName);

        // CollectMultipleMetricsSpark provider that creates an initializes a custom
        // collector
        CollectMultipleMetricsSpark.SparkCollectorProvider customProvider =
            new CollectMultipleMetricsSpark.SparkCollectorProvider() {
                @Override
                public MetricsCollectorSpark<? extends MetricsArgumentCollection> createCollector(
                    final String outputBaseName,
                    final Set<MetricAccumulationLevel> metricAccumulationLevel,
                    final List<Header> defaultHeaders,
                    final SAMFileHeader samHeader)
                {
                    return testCollector;
                }
        };

        // Manually create a tool and programmatically set the custome collector as the one
        // to run
        CollectMultipleMetricsSpark multipleCollectorTool = new CollectMultipleMetricsSpark();
        multipleCollectorTool.setCollectorsToRun(Collections.singletonList(customProvider));
        multipleCollectorTool.instanceMain(args.getArgsArray());

        Assert.assertEquals(testCollector.count, 52L);
    }

}
