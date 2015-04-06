package org.broadinstitute.hellbender.tools.picard.illumina.metrics;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.util.*;

public class CollectIlluminaBasecallingMetricsTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/illumina/metrics/CollectIlluminaBasecallingMetrics");

    private File rootTestDir;

    @BeforeTest
    private void setUp() throws Exception {
        rootTestDir = File.createTempFile("cibm.", ".tmp");
        Assert.assertTrue(rootTestDir.delete());
        Assert.assertTrue(rootTestDir.mkdir());
        for (final File source : TEST_DATA_DIR.listFiles()) {
            if (source.isDirectory() && !source.isHidden()) {
                IOUtil.copyDirectoryTree(source, new File(rootTestDir.getPath(), source.getName()));
            }
        }
    }

    @AfterTest
    private void tearDown() {
        IOUtil.deleteDirectoryTree(rootTestDir);
    }

    @Test
    public void testIndexedRunLane1() throws Exception {
        final MetricsFile<IlluminaBasecallingMetrics, Integer> metricsFile = runIt(1, "25T8B25T","25T8B25T/Data/Intensities/BaseCalls", true);
        final IlluminaBasecallingMetrics metric1 = metricsFile.getMetrics().get(0);
        Assert.assertEquals(metric1.LANE, "1");
        Assert.assertEquals(metric1.MOLECULAR_BARCODE_SEQUENCE_1, "AACAATGG");
        Assert.assertEquals(metric1.MOLECULAR_BARCODE_NAME, "tagged_117");
        Assert.assertEquals(metric1.MEAN_CLUSTERS_PER_TILE, 2.0);
        Assert.assertEquals(metric1.SD_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(metric1.MEAN_PF_CLUSTERS_PER_TILE, 2.0);
        Assert.assertEquals(metric1.SD_PF_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(metric1.MEAN_PCT_PF_CLUSTERS_PER_TILE, 100.0);
        Assert.assertEquals(metric1.SD_PCT_PF_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(metric1.TOTAL_CLUSTERS, 2);
        Assert.assertEquals(metric1.TOTAL_CLUSTERS * 50, metric1.TOTAL_BASES);
        Assert.assertEquals(metric1.PF_BASES, metric1.TOTAL_BASES * metric1.PF_CLUSTERS / metric1.TOTAL_CLUSTERS);

        final IlluminaBasecallingMetrics metric2 = metricsFile.getMetrics().get(1);
        Assert.assertEquals(metric2.LANE, "1");
        Assert.assertEquals(metric2.MOLECULAR_BARCODE_SEQUENCE_1, "AACGCATT");
        Assert.assertEquals(metric2.MOLECULAR_BARCODE_NAME, "tagged_741");
        Assert.assertEquals(metric2.MEAN_CLUSTERS_PER_TILE, 3.0);
        Assert.assertEquals(metric2.SD_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(metric2.MEAN_PF_CLUSTERS_PER_TILE, 2.0);
        Assert.assertEquals(metric2.SD_PF_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(metric2.MEAN_PCT_PF_CLUSTERS_PER_TILE, 66.67);
        Assert.assertEquals(metric2.SD_PCT_PF_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(metric2.TOTAL_CLUSTERS, 3);
        Assert.assertEquals(metric2.TOTAL_CLUSTERS * 50, metric2.TOTAL_BASES);
        Assert.assertEquals(metric2.PF_BASES, metric2.TOTAL_BASES * metric2.PF_CLUSTERS / metric2.TOTAL_CLUSTERS);

        final IlluminaBasecallingMetrics metric3 = metricsFile.getMetrics().get(2);
        Assert.assertEquals(metric3.LANE, "1");
        Assert.assertEquals(metric3.MOLECULAR_BARCODE_SEQUENCE_1, "ACAGGTAT");
        Assert.assertEquals(metric3.MOLECULAR_BARCODE_NAME, "tagged_375");
        Assert.assertEquals(metric3.MEAN_CLUSTERS_PER_TILE, 1.0);
        Assert.assertEquals(metric3.SD_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(metric3.MEAN_PF_CLUSTERS_PER_TILE, 1.0);
        Assert.assertEquals(metric3.SD_PF_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(metric3.MEAN_PCT_PF_CLUSTERS_PER_TILE, 100.0);
        Assert.assertEquals(metric3.SD_PCT_PF_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(metric3.TOTAL_CLUSTERS, 1);
        Assert.assertEquals(metric3.TOTAL_CLUSTERS * 50, metric3.TOTAL_BASES);
        Assert.assertEquals(metric3.PF_BASES, metric3.TOTAL_BASES * metric3.PF_CLUSTERS / metric3.TOTAL_CLUSTERS);

        final IlluminaBasecallingMetrics metric4 = metricsFile.getMetrics().get(3);
        Assert.assertEquals(metric4.LANE, "1");
        Assert.assertEquals(metric4.MOLECULAR_BARCODE_SEQUENCE_1, "ACTAAGAC");
        Assert.assertEquals(metric4.MOLECULAR_BARCODE_NAME, "tagged_630");
        Assert.assertEquals(metric4.MEAN_CLUSTERS_PER_TILE, 2.0);
        Assert.assertEquals(metric4.SD_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(metric4.MEAN_PF_CLUSTERS_PER_TILE, 1.0);
        Assert.assertEquals(metric4.SD_PF_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(metric4.MEAN_PCT_PF_CLUSTERS_PER_TILE, 50.00);
        Assert.assertEquals(metric4.SD_PCT_PF_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(metric4.TOTAL_CLUSTERS, 2);
        Assert.assertEquals(metric4.TOTAL_CLUSTERS * 50, metric4.TOTAL_BASES);
        Assert.assertEquals(metric4.PF_BASES, metric4.TOTAL_BASES * metric4.PF_CLUSTERS / metric4.TOTAL_CLUSTERS);

        final IlluminaBasecallingMetrics metric5 = metricsFile.getMetrics().get(4);
        Assert.assertEquals(metric5.LANE, "1");
        Assert.assertEquals(metric5.MOLECULAR_BARCODE_SEQUENCE_1, "AGCATGGA");
        Assert.assertEquals(metric5.MOLECULAR_BARCODE_NAME, "tagged_908");
        Assert.assertEquals(metric5.MEAN_CLUSTERS_PER_TILE, 1.0);
        Assert.assertEquals(metric5.SD_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(metric5.MEAN_PF_CLUSTERS_PER_TILE, 1.0);
        Assert.assertEquals(metric5.SD_PF_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(metric5.MEAN_PCT_PF_CLUSTERS_PER_TILE, 100.0);
        Assert.assertEquals(metric5.SD_PCT_PF_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(metric5.TOTAL_CLUSTERS, 1);
        Assert.assertEquals(metric5.TOTAL_CLUSTERS * 50, metric5.TOTAL_BASES);
        Assert.assertEquals(metric5.PF_BASES, metric5.TOTAL_BASES * metric5.PF_CLUSTERS / metric5.TOTAL_CLUSTERS);

        final IlluminaBasecallingMetrics laneMetric = metricsFile.getMetrics().get(34);
        Assert.assertEquals(laneMetric.LANE, "1");
        Assert.assertEquals(laneMetric.MOLECULAR_BARCODE_SEQUENCE_1, null);
        Assert.assertEquals(laneMetric.MOLECULAR_BARCODE_NAME, null);
        Assert.assertEquals(laneMetric.MEAN_CLUSTERS_PER_TILE, 60.0);
        Assert.assertEquals(laneMetric.SD_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(laneMetric.MEAN_PF_CLUSTERS_PER_TILE, 50.0);
        Assert.assertEquals(laneMetric.SD_PF_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(laneMetric.MEAN_PCT_PF_CLUSTERS_PER_TILE, 83.33);
        Assert.assertEquals(laneMetric.SD_PCT_PF_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(laneMetric.TOTAL_CLUSTERS, 60);
        Assert.assertEquals(laneMetric.TOTAL_CLUSTERS * 50, laneMetric.TOTAL_BASES);
        Assert.assertEquals(laneMetric.PF_BASES, laneMetric.TOTAL_BASES * laneMetric.PF_CLUSTERS / laneMetric.TOTAL_CLUSTERS);
    }

    @Test
    public void testNonIndexedRunLane1() throws Exception {
        final MetricsFile<IlluminaBasecallingMetrics, Integer> metricsFile = runIt(1, "125T125T","125T125T/Data/Intensities/BaseCalls",false);
        final IlluminaBasecallingMetrics laneMetric = metricsFile.getMetrics().get(0);

        Assert.assertEquals(laneMetric.LANE, "1");
        Assert.assertEquals(laneMetric.MOLECULAR_BARCODE_SEQUENCE_1, null);
        Assert.assertEquals(laneMetric.MOLECULAR_BARCODE_NAME, null);
        Assert.assertEquals(laneMetric.MEAN_CLUSTERS_PER_TILE, 2000.0);
        Assert.assertEquals(laneMetric.SD_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(laneMetric.MEAN_PF_CLUSTERS_PER_TILE, 1863.0);
        Assert.assertEquals(laneMetric.SD_PF_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(laneMetric.MEAN_PCT_PF_CLUSTERS_PER_TILE, 93.15);
        Assert.assertEquals(laneMetric.SD_PCT_PF_CLUSTERS_PER_TILE, 0.0);
        Assert.assertEquals(laneMetric.TOTAL_BASES, 500000);
        Assert.assertEquals(laneMetric.TOTAL_READS, 4000);
        Assert.assertEquals(laneMetric.PF_BASES, 465750);
        Assert.assertEquals(laneMetric.PF_READS, 3726);


        Assert.assertEquals(metricsFile.getMetrics().size(), 1);
    }

    private MetricsFile<IlluminaBasecallingMetrics, Integer> runIt(final int lane, final String readStructure, final String basecallsDirName, final boolean isIndexed) throws Exception {
        final File metricsFile = File.createTempFile("cibm.", ".metrics");
        metricsFile.deleteOnExit();

        File basecallsDir = new File(rootTestDir.getPath(),basecallsDirName);

        List<String> argsList = new ArrayList<String>();
        argsList.add("--BASECALLS_DIR"); argsList.add(basecallsDir.getPath());
        argsList.add("--LANE"); argsList.add(Integer.toString(lane));
        argsList.add("--OUTPUT"); argsList.add(metricsFile.getPath());

        if (readStructure != null) {
            argsList.add("--READ_STRUCTURE"); argsList.add(readStructure);
        }
        if (isIndexed) {
            argsList.add("--INPUT"); argsList.add(new File(basecallsDir.getPath(),"barcodeData." + lane).getPath());
        }

        runCommandLine(argsList);

        final MetricsFile<IlluminaBasecallingMetrics,Integer> retval = new MetricsFile<IlluminaBasecallingMetrics,Integer>();
        retval.read(new FileReader(metricsFile));
        return retval;
    }
}
