package org.broadinstitute.hellbender.tools.spark.pathseq.loggers;

import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;

public class PSFilterFileLoggerTest  extends GATKBaseTest {
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
        filterLogger.close();

        final File expectedOutputFile = getTestFile("expected.filter.metrics");

        Assert.assertTrue(MetricsFile.areMetricsEqual(metricsOutputFile, expectedOutputFile));
    }
}