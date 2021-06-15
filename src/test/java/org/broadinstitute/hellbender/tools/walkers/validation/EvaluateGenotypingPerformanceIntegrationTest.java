package org.broadinstitute.hellbender.tools.walkers.validation;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
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
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth.vcf", new File(TEST_DIR + "expected_correlations_basic.tsv"), new File(TEST_DIR + "expected_accuracy_basic.tsv"),
                        Collections.emptyList()},
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth_no_monomorphic.vcf", new File(TEST_DIR + "expected_correlations_basic.tsv"), new File(TEST_DIR + "expected_accuracy_basic.tsv"),
                        Collections.emptyList()}, //sites not included in truth should be treated as truth being hom ref for all samples
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth.vcf", new File(TEST_DIR + "expected_correlations_dosage.tsv"), new File(TEST_DIR + "expected_accuracy_basic.tsv"),
                        Arrays.asList("--dosage-field", "DS")},
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth.vcf", new File(TEST_DIR + "expected_correlations_basic.tsv"), new File(TEST_DIR + "expected_accuracy_min_af.tsv"),
                        Arrays.asList("--min-af-for-accuracy-metrics", "0.2")},
                {TEST_DIR + "eval_X.vcf", TEST_DIR + "truth_X.vcf", new File(TEST_DIR + "expected_correlations_X.tsv"), new File(TEST_DIR + "expected_accuracy_X.tsv"),
                Arrays.asList("--allow-differing-ploidies")},
                {TEST_DIR + "eval_extra_sample.vcf", TEST_DIR + "truth.vcf", new File(TEST_DIR + "expected_correlations_basic.tsv"), new File(TEST_DIR + "expected_accuracy_basic.tsv"),
                        Arrays.asList("--af-annotations", "SAMPLE3:AF3")}, //extra sample only in eval should be ignored
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth.vcf", new File(TEST_DIR + "expected_correlations_af_resource.tsv"), new File(TEST_DIR + "expected_accuracy_basic.tsv"),
                        Arrays.asList("--resource", TEST_DIR + "af_resource.vcf")},
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth_new_sample_names.vcf", new File(TEST_DIR + "expected_correlations_basic.tsv"), new File(TEST_DIR + "expected_accuracy_basic.tsv"),
                        Arrays.asList("--sample-map", "SAMPLE1:SAMPLE3", "--sample-map", "SAMPLE2:SAMPLE4")},

        };
    }

    @Test(dataProvider = "integrationTestDataProvider")
    void integrationTest(final String eval, final String truth, final File expectedCorrelation, final File expectedAccuracy, final List<String> additionalArgs) throws IOException {
        final File outputCorrelation = createTempFile("correlation", ".tsv");
        final File outputAccuracy = createTempFile("accuracy", ".tsv");

        runCommandLine(buildArgs(eval, truth, outputCorrelation.getAbsolutePath(), outputAccuracy.getAbsolutePath(), additionalArgs));

        IntegrationTestSpec.assertEqualTextFiles(outputCorrelation, expectedCorrelation, "#");
        IntegrationTestSpec.assertEqualTextFiles(outputAccuracy, expectedAccuracy, "#");
    }

    @Test(expectedExceptions = GATKException.class)
    void differingPloidyTest() {
        runCommandLine(buildArgs(TEST_DIR + "eval_X.vcf", TEST_DIR + "truth_X.vcf",
                createTempFile("correlation", ".tsv").getAbsolutePath(), createTempFile("accuracy", ".tsv").getAbsolutePath(),
                Collections.emptyList()));
    }

    List<String> buildArgs(final String evalPath, final String truthPath, final String outputCorrelationPath, final String outputAccuracyPath, final List<String> additionalArgs) {
        final List<String> args = new ArrayList<>(Arrays.asList(
                "--first-bin-right-edge", "0.1",
                "-nbins", "5",
                "--af-annotations", "SAMPLE1:AF1",
                "--af-annotations", "SAMPLE2:AF2",
                "--eval", evalPath,
                "--truth", truthPath,
                "--O", outputCorrelationPath,
                "--OA", outputAccuracyPath
        ));
        args.addAll(additionalArgs);

        return args;
    }
}
