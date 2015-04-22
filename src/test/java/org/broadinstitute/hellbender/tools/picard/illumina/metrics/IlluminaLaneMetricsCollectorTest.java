package org.broadinstitute.hellbender.tools.picard.illumina.metrics;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.picard.illumina.parser.ReadStructure;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

/** @author mccowan */
public final class IlluminaLaneMetricsCollectorTest extends CommandLineProgramTest {
    final static File TEST_DIRECTORY = new File(getTestDataDir(), "picard/illumina/metrics/IlluminaLaneMetricsCollector");

    private static File buildOutputFile(final File directory, final String prefix, final String extension) {
        return new File(directory, String.format("%s.%s", prefix, extension));
    }

    @Test(dataProvider = "testCollectIlluminaLaneMetrics")
    public void testCollectIlluminaLaneMetrics(final String testRun, final ReadStructure readStructure) throws Exception {
        final File runDirectory = new File(TEST_DIRECTORY, testRun);
        final CollectIlluminaLaneMetrics clp = new CollectIlluminaLaneMetrics();
        clp.OUTPUT_DIRECTORY = IOUtil.createTempDir("illuminaLaneMetricsCollectorTest", null);
        clp.RUN_DIRECTORY = runDirectory;
        clp.OUTPUT_PREFIX = testRun;
        clp.READ_STRUCTURE = readStructure;
        clp.doWork();

        final File phasingMetricsPhile = buildOutputFile(clp.OUTPUT_DIRECTORY, clp.OUTPUT_PREFIX, IlluminaPhasingMetrics.getExtension());
        final File canonicalPhasingPhile = buildOutputFile(runDirectory, testRun, IlluminaPhasingMetrics.getExtension());
        Assert.assertTrue(MetricsFile.areMetricsEqual(canonicalPhasingPhile, phasingMetricsPhile));

        final File laneMetricsFile = buildOutputFile(clp.OUTPUT_DIRECTORY, clp.OUTPUT_PREFIX, IlluminaLaneMetrics.getExtension());
        final File canonicalLaneFile = buildOutputFile(runDirectory, testRun, IlluminaLaneMetrics.getExtension());
        Assert.assertTrue(MetricsFile.areMetricsEqual(canonicalLaneFile, laneMetricsFile));

        IOUtil.deleteDirectoryTree(clp.OUTPUT_DIRECTORY);
    }

    @DataProvider(name = "testCollectIlluminaLaneMetrics")
    public Object[][] testCollectIlluminaLaneMetricsDataProvider() {
        return new Object[][] {
                {"A7LE0", new ReadStructure("25T8B8B25T")},
                {"C2MFAACXX", new ReadStructure("95T101T")},
                {"H7BATADXX", new ReadStructure("76T8B76T")},
                {"H7H7RADXX", new ReadStructure("101T8B8B101T")},
                {"A67HY", new ReadStructure("8B8B")}
        };
    }
}
