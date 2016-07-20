package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.ValidationStringency;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.metrics.InsertSizeMetricsArgumentCollection;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.tools.spark.pipelines.metrics.InsertSizeMetricsCollectorSpark;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;

/**
 * This test uses the same input and output as CollectInsertSizeMetricsSparkIntegrationTest,
 * except that it bypasses CollectInsertSizeMetricsSpark and uses InsertSizeMetricsCollectorSpark
 * directly in order to force the input RDD to be split across two partitions. This ensures that
 * the MultiLevelCollector combine and combineUnit methods are executed during the test.
 */
public class InsertSizeMetricsCollectorSparkUnitTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/analysis/CollectInsertSizeMetrics");

    public String getTestedClassName() {
        return InsertSizeMetricsCollectorSpark.class.getSimpleName();
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

    @Test(dataProvider="metricsfiles", groups="spark")
    public void test(
            final String fileName,
            final String referenceName,
            final boolean allLevels,
            final String expectedResultsFile) throws IOException {

        final String inputPath = new File(TEST_DATA_DIR, fileName).getAbsolutePath();
        final String referencePath = referenceName != null ? new File(referenceName).getAbsolutePath() : null;

        final File outfile = BaseTest.createTempFile("test", ".insert_size_metrics");

        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        ReadsSparkSource readSource = new ReadsSparkSource(ctx, ValidationStringency.DEFAULT_STRINGENCY);

        SAMFileHeader samHeader = readSource.getHeader(inputPath, referencePath, null);
        JavaRDD<GATKRead> rddParallelReads = readSource.getParallelReads(inputPath, referencePath);

        InsertSizeMetricsArgumentCollection isArgs = new InsertSizeMetricsArgumentCollection();
        isArgs.output = outfile.getAbsolutePath();
        isArgs.useEnd = InsertSizeMetricsArgumentCollection.EndToUse.SECOND;
        if (allLevels) {
            isArgs.metricAccumulationLevel = new HashSet<>();
            isArgs.metricAccumulationLevel.add(MetricAccumulationLevel.ALL_READS);
            isArgs.metricAccumulationLevel.add(MetricAccumulationLevel.SAMPLE);
            isArgs.metricAccumulationLevel.add(MetricAccumulationLevel.LIBRARY);
            isArgs.metricAccumulationLevel.add(MetricAccumulationLevel.READ_GROUP);
        }

        InsertSizeMetricsCollectorSpark isSpark = new InsertSizeMetricsCollectorSpark();
        isSpark.initialize(isArgs, samHeader, null);
        ReadFilter rf = isSpark.getReadFilter(samHeader);

        // Force the input RDD to be split into two partitions to ensure that the
        // reduce/combiners run
        rddParallelReads = rddParallelReads.repartition(2);
        isSpark.collectMetrics(rddParallelReads.filter(r -> rf.test(r)), samHeader);

        isSpark.saveMetrics(fileName, null);

        IntegrationTestSpec.assertEqualTextFiles(
                outfile,
                new File(TEST_DATA_DIR, expectedResultsFile),
                "#"
        );
    }

}
