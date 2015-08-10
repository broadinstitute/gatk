package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
    import org.broadinstitute.hellbender.CommandLineProgramTest;
    import org.broadinstitute.hellbender.tools.IntegrationTestSpec;
    import org.testng.Assert;
    import org.testng.annotations.Test;

    import java.io.File;
    import java.io.FileReader;
    import java.io.IOException;

public final class CollectBaseDistributionByCycleIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/analysis/CollectBaseDistributionByCycle");

    @Test(groups = {"R"})
    public void test() throws IOException {
        final File input = new File(TEST_DATA_DIR, "first5000a.bam");
        final File expectedFile = new File(TEST_DATA_DIR, "CollectBaseDistributionByCycle.txt");
        final File outfile = createTempFile("testCollectBaseDistributionByCycle", ".metrics");
        final File pdf = createTempFile("testCollectBaseDistributionByCycle", ".pdf");
        outfile.deleteOnExit();
        pdf.deleteOnExit();
        final String[] args = new String[]{
                "--INPUT", input.getAbsolutePath(),
                "--OUTPUT", outfile.getAbsolutePath(),
                "--CHART", pdf.getAbsolutePath(),
                "--PRODUCE_PLOT", "true",
        };
        runCommandLine(args);

        try (final FileReader actualReader = new FileReader(outfile);) {
            final MetricsFile<?,Integer> output = new MetricsFile<>();
            output.read(actualReader);
            Assert.assertEquals(output.getMetrics().size(), 202);
        }
        Assert.assertTrue(pdf.exists());
        Assert.assertTrue(pdf.length() > 0);
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }
}