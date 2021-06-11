package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class EvaluateGenotypingPerformanceIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_DIR = getTestDataDir() + "/EvaluateGenotypingPerformance/";

    @DataProvider(name="integrationTestDataProvider")
    Object[][] integrationTestDataProvider() {
        return new Object[][]{
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth.vcf", new File(TEST_DIR + "expected_correlations_basic.tsv"), new File(TEST_DIR + "expected_accuracy_basic.tsv"), Collections.emptyList()},
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth.vcf", new File(TEST_DIR + "expected_correlations_dosage.tsv"), new File(TEST_DIR + "expected_accuracy_basic.tsv"),
                        Arrays.asList("--dosage-field", "DS")},
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth.vcf", new File(TEST_DIR + "expected_correlations_basic.tsv"), new File(TEST_DIR + "expected_accuracy_min_af.tsv"),
                        Arrays.asList("--min-af-for-accuracy-metrics", "0.2")}
        };
    }

    @Test(dataProvider = "integrationTestDataProvider")
    void integrationTest(final String eval, final String truth, final File expectedCorrelation, final File expectedAccuracy, final List<String> additionalArgs) throws IOException {
        final File outputCorrelation = createTempFile("correlation", ".tsv");
        final File outputAccuracy = createTempFile("accuracy", ".tsv");
        final List<String> args = new ArrayList<>(Arrays.asList(
                "--first-bin-right-edge", "0.1",
                "-nbins", "5",
                "--af-annotations", "SAMPLE1:AF1",
                "--af-annotations", "SAMPLE2:AF2",
                "--eval", eval,
                "--truth", truth,
                "--O", outputCorrelation.getAbsolutePath(),
                "--OA", outputAccuracy.getAbsolutePath()
        ));
        args.addAll(additionalArgs);

        runCommandLine(args);
        IntegrationTestSpec.assertEqualTextFiles(outputCorrelation, expectedCorrelation, "#");
        IntegrationTestSpec.assertEqualTextFiles(outputAccuracy, expectedAccuracy, "#");
    }
}
