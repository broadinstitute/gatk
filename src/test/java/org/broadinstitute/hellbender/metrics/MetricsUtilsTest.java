package org.broadinstitute.hellbender.metrics;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.hadoop.hdfs.MiniDFSCluster;
import org.broadinstitute.hellbender.utils.gcs.GoogleStorageUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.MiniClusterUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

public class MetricsUtilsTest extends GATKBaseTest {
    private MiniDFSCluster cluster;
    private String hdfsWorkingDir;

    @BeforeClass(alwaysRun = true)
    private void setupMiniCluster() throws IOException {
        cluster = MiniClusterUtils.getMiniCluster();
        hdfsWorkingDir = MiniClusterUtils.getWorkingDir(cluster).toString();
    }

    @AfterClass(alwaysRun = true)
    private void shutdownMiniCluster() {
       MiniClusterUtils.stopCluster(cluster);
    }

    @DataProvider(name = "metricsPaths")
    public Object[][] getMetricsPaths(){
        return new Object[][]{
                {"metrics"},
                {getGCPTestStaging()},
                {hdfsWorkingDir}
        };
    }

    public static class TestMetric extends MetricBase {
        public Integer value1 = 0;
        public Integer value2 = 0;
    }

    @Test(dataProvider = "metricsPaths", groups = "bucket")
    public void testSaveMetrics(String destinationPrefix) throws IOException {
        final String outputPath = GoogleStorageUtils.getTempFilePath(destinationPrefix, ".txt");
        TestMetric testMetric = new TestMetric();
        testMetric.value1 = 10;
        testMetric.value2 = 5;

        final MetricsFile<TestMetric, ?> metrics = new MetricsFile<>();
        metrics.addMetric(testMetric);
        MetricsUtils.saveMetrics(metrics, outputPath);
        Assert.assertTrue(Files.exists(IOUtils.getPath(outputPath)));
        File localCopy = copyFileToLocalTmpFile(outputPath);

        final File expectedMetrics = createTempFile("expectedMetrics", ".txt");
        metrics.write(expectedMetrics);

        Assert.assertTrue(MetricsFile.areMetricsEqual(localCopy, expectedMetrics));
    }

    private File copyFileToLocalTmpFile(String outputPath) throws IOException {
        File localCopy = createTempFile("local_metrics_copy",".txt");
        Files.copy(IOUtils.getPath(outputPath), IOUtils.getPath(localCopy.getAbsolutePath()), StandardCopyOption.REPLACE_EXISTING);
        return localCopy;
    }
}
