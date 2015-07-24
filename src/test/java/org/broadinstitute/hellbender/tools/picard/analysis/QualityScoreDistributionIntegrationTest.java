package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class QualityScoreDistributionIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/analysis/QualityScoreDistribution");

    @Test
    public void test() throws IOException {
        final File input = new File(TEST_DATA_DIR, "first5000a.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "qualscoredist.txt");
        final File outfile = BaseTest.createTempFile("testQualityScoreDistribution", ".metrics");
        final File pdf = BaseTest.createTempFile("testQualityScoreDistribution", ".pdf");
        final String[] args = new String[]{
                "--INPUT", input.getAbsolutePath(),
                "--OUTPUT", outfile.getAbsolutePath(),
                "--CHART", pdf.getAbsolutePath(),
                "--PRODUCE_PLOT", "true",
        };
        runCommandLine(args);

        try (final FileReader actualReader = new FileReader(outfile);) {
            final MetricsFile<?, Integer> output = new MetricsFile<>();
            output.read(actualReader);
            Assert.assertEquals(output.getAllHistograms().size(), 1);
            Assert.assertEquals(output.getHistogram().size(), 41);
        }
        Assert.assertTrue(pdf.exists());
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }
}
