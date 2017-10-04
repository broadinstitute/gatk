package org.broadinstitute.hellbender.tools.spark.pathseq.loggers;

import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.metrics.MetricsUtils;
import org.broadinstitute.hellbender.tools.spark.pathseq.PSScorer;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class PSScoreFileLoggerTest extends BaseTest {

    @Test
    public void test() {
        final MetricsFile<PSScoreMetrics, Long> metricsFile = new MetricsFile<>();
        final File metricsOutputFile = createTempFile("metrics",".txt");
        final PSScoreFileLogger scoreLogger = new PSScoreFileLogger(metricsFile, metricsOutputFile.getAbsolutePath());

        final int numMappedReads = 21;
        final int numUnmappedReads = 7;
        final List<GATKRead> reads = new ArrayList<>(numMappedReads + numUnmappedReads);
        for (int i = 0; i < numMappedReads; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead("101M");
            read.setAttribute(PSScorer.HITS_TAG, "42");
            reads.add(read);
        }
        for (int i = 0; i < numUnmappedReads; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead("101M");
            reads.add(read);
        }

        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final JavaRDD<GATKRead> readRdd = ctx.parallelize(reads);
        scoreLogger.logReadCounts(readRdd);
        scoreLogger.writeFile();

        final MetricsFile<PSScoreMetrics, Long> expectedMetricsFile = new MetricsFile<>();
        final PSScoreMetrics expectedScoreMetrics = new PSScoreMetrics();
        expectedScoreMetrics.UNMAPPED_READS = (long) numUnmappedReads;
        expectedScoreMetrics.MAPPED_READS = (long) numMappedReads;
        expectedMetricsFile.addMetric(expectedScoreMetrics);
        final File expectedOutputFile = createTempFile("expected_metrics",".txt");
        MetricsUtils.saveMetrics(expectedMetricsFile, expectedOutputFile.getAbsolutePath());

        Assert.assertTrue(MetricsFile.areMetricsEqual(metricsOutputFile, expectedOutputFile));
    }

}