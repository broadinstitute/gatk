package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public final class CollectQualityYieldMetricsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/analysis/CollectQualityYieldMetrics");

    @Test
    public void test() throws IOException {
        final File input = new File(TEST_DATA_DIR, "collect_quality_yield_metrics.sam");
        final File expectedFile = new File(TEST_DATA_DIR, "collect_quality_yield_metrics.txt");   //file created using picard 1.130
        final File outfile = BaseTest.createTempFile("testCollectQualityYield", ".metrics");
        final String[] args = new String[]{
                "--INPUT", input.getAbsolutePath(),
                "--OUTPUT", outfile.getAbsolutePath(),
        };
        runCommandLine(args);

        try (final FileReader actualReader = new FileReader(outfile);) {
            final MetricsFile<?,Integer> output = new MetricsFile<>();
            output.read(actualReader);
            Assert.assertEquals(output.getMetrics().size(), 1);
            Assert.assertEquals(output.getAllHistograms().size(), 0);
        }
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }
}
