package org.broadinstitute.hellbender.tools.spark.pathseq.loggers;

import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.metrics.MetricsUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;

public class PSFilterFileLoggerTest  extends BaseTest {

    @Test
    public void test() {
        final MetricsFile<PSFilterMetrics, Long> metricsFile = new MetricsFile<>();
        final File metricsOutputFile = createTempFile("metrics", ".txt");
        final PSFilterFileLogger filterLogger = new PSFilterFileLogger(metricsFile, metricsOutputFile.getAbsolutePath());

        final int numPrimaryReads = 1000;
        final int numReadsAfterPrealign = 950;
        final int numReadsAfterQuality = 875;
        final int numReadsAfterHost = 775;
        final int numReadsAfterDuplicate = 760;
        final int numFinalPaired = 300;
        final GATKRead read = ArtificialReadUtils.createArtificialRead("101M");

        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        filterLogger.logPrimaryReads(ctx.parallelize(Collections.nCopies(numPrimaryReads, read)));
        filterLogger.logReadsAfterPrealignedHostFilter(ctx.parallelize(Collections.nCopies(numReadsAfterPrealign, read)));
        filterLogger.logReadsAfterQualityFilter(ctx.parallelize(Collections.nCopies(numReadsAfterQuality, read)));
        filterLogger.logReadsAfterHostFilter(ctx.parallelize(Collections.nCopies(numReadsAfterHost, read)));
        filterLogger.logReadsAfterDeduplication(ctx.parallelize(Collections.nCopies(numReadsAfterDuplicate, read)));
        filterLogger.logFinalPairedReads(ctx.parallelize(Collections.nCopies(numFinalPaired, read)));
        filterLogger.writeFile();

        final MetricsFile<PSFilterMetrics, Long> expectedMetricsFile = new MetricsFile<>();
        final PSFilterMetrics expectedFilterMetrics = new PSFilterMetrics();
        expectedFilterMetrics.PRIMARY_READS = (long) numPrimaryReads;
        expectedFilterMetrics.READS_AFTER_PREALIGNED_HOST_FILTER = (long) numReadsAfterPrealign;
        expectedFilterMetrics.READS_AFTER_QUALITY_AND_COMPLEXITY_FILTER = (long) numReadsAfterQuality;
        expectedFilterMetrics.READS_AFTER_HOST_FILTER = (long) numReadsAfterHost;
        expectedFilterMetrics.READS_AFTER_DEDUPLICATION = (long) numReadsAfterDuplicate;
        expectedFilterMetrics.FINAL_PAIRED_READS = (long) numFinalPaired;

        //Derived metrics
        expectedFilterMetrics.FINAL_TOTAL_READS = expectedFilterMetrics.READS_AFTER_DEDUPLICATION;
        expectedFilterMetrics.FINAL_UNPAIRED_READS = expectedFilterMetrics.FINAL_TOTAL_READS - expectedFilterMetrics.FINAL_PAIRED_READS;
        expectedFilterMetrics.LOW_QUALITY_OR_LOW_COMPLEXITY_READS_FILTERED = expectedFilterMetrics.READS_AFTER_PREALIGNED_HOST_FILTER - expectedFilterMetrics.READS_AFTER_QUALITY_AND_COMPLEXITY_FILTER;
        expectedFilterMetrics.HOST_READS_FILTERED = expectedFilterMetrics.PRIMARY_READS - expectedFilterMetrics.READS_AFTER_PREALIGNED_HOST_FILTER + expectedFilterMetrics.READS_AFTER_QUALITY_AND_COMPLEXITY_FILTER - expectedFilterMetrics.READS_AFTER_HOST_FILTER;
        expectedFilterMetrics.DUPLICATE_READS_FILTERED = expectedFilterMetrics.READS_AFTER_HOST_FILTER - expectedFilterMetrics.READS_AFTER_DEDUPLICATION;

        expectedMetricsFile.addMetric(expectedFilterMetrics);
        final File expectedOutputFile = createTempFile("expected_metrics", ".txt");
        MetricsUtils.saveMetrics(expectedMetricsFile, expectedOutputFile.getAbsolutePath());

        Assert.assertTrue(MetricsFile.areMetricsEqual(metricsOutputFile, expectedOutputFile));
    }
}