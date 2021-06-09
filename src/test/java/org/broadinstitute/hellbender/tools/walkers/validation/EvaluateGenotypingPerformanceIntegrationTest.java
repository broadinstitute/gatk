package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class EvaluateGenotypingPerformanceIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_DIR = getTestDataDir() + "/EvaluateGenotypingPerformance/";

    @Test
    void basicTest() {
        final File outputCorrelation = createTempFile("correlation", ".tsv");
        final File outputAccuracy = createTempFile("accuracy", ".tsv");
        final String[] args = {
                "--first-bin-right-edge", "0.1",
                "-nbins", "5",
                "--af-annotations", "SAMPLE1:AF1",
                "--af-annotations", "SAMPLE2:AF2",
                "--eval", TEST_DIR + "eval.vcf",
                "--truth", TEST_DIR + "truth.vcf",
                "--O", outputCorrelation.getAbsolutePath(),
                "--OA", outputAccuracy.getAbsolutePath()
        };

        runCommandLine(args);

        assertMetricsAreEqual(outputCorrelation, new File(TEST_DIR + "expected_correlations_basic.tsv"));
    }

    void assertMetricsAreEqual(final File actualMetricsFile, final File expectedMetricsFile) {
        final MetricsFile<?, ?> actualMetrics = new MetricsFile<>();
        try (final FileReader reader = new FileReader(actualMetricsFile)) {
            actualMetrics.read(reader);
        } catch (final IOException ex) {
            throw new GATKException("Error reading metrics file " + actualMetricsFile, ex);
        }

        final MetricsFile<?, ?> expectedMetrics = new MetricsFile<>();
        try (final FileReader reader = new FileReader(expectedMetricsFile)) {
            expectedMetrics.read(reader);
        } catch (final IOException ex) {
            throw new GATKException("Error reading metrics file " + expectedMetrics, ex);
        }

        Assert.assertEquals(actualMetrics.getMetrics(), expectedMetrics.getMetrics());
    }
}
