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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class EvaluateGenotypingPerformanceIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_DIR = getTestDataDir() + "/EvaluateGenotypingPerformance/";

    @DataProvider(name="integrationTestDataProvider")
    Object[][] integrationTestDataProvider() {
        final File expectedAccuracyAFBasic = new File(TEST_DIR + "expected_accuracy_af_basic.tsv");
        return new Object[][]{
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth.vcf", new File(TEST_DIR + "expected_correlations_basic.tsv"), new File(TEST_DIR + "expected_accuracy_basic.tsv"),
                        expectedAccuracyAFBasic, Collections.emptyList(), Collections.emptyMap()},
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth_no_monomorphic.vcf", new File(TEST_DIR + "expected_correlations_basic.tsv"), new File(TEST_DIR + "expected_accuracy_basic.tsv"),
                        expectedAccuracyAFBasic, Collections.emptyList(), Collections.emptyMap()}, //sites not included in truth should be treated as truth being hom ref for all samples
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth.vcf", new File(TEST_DIR + "expected_correlations_dosage.tsv"), new File(TEST_DIR + "expected_accuracy_basic.tsv"),
                        expectedAccuracyAFBasic, Arrays.asList("--dosage-field", "DS"), Collections.emptyMap()},
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth.vcf", new File(TEST_DIR + "expected_correlations_basic.tsv"), new File(TEST_DIR + "expected_accuracy_min_af.tsv"),
                        expectedAccuracyAFBasic, Arrays.asList("--min-af-for-accuracy-metrics", "0.2"), Collections.emptyMap()}, //min af does not apply to af-stratified accuracy
                {TEST_DIR + "eval_X.vcf", TEST_DIR + "truth_X.vcf", new File(TEST_DIR + "expected_correlations_X.tsv"), new File(TEST_DIR + "expected_accuracy_X.tsv"),
                        new File(TEST_DIR + "expected_accuracy_af_X.tsv"), Arrays.asList("--allow-differing-ploidies"), Collections.emptyMap()},
                {TEST_DIR + "eval_extra_sample.vcf", TEST_DIR + "truth.vcf", new File(TEST_DIR + "expected_correlations_basic.tsv"), new File(TEST_DIR + "expected_accuracy_basic.tsv"),
                        expectedAccuracyAFBasic, Arrays.asList("--af-annotations", "SAMPLE3:AF3"), Collections.emptyMap()}, //extra sample only in eval should be ignored
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth.vcf", new File(TEST_DIR + "expected_correlations_af_resource.tsv"), new File(TEST_DIR + "expected_accuracy_basic.tsv"),
                        new File(TEST_DIR + "expected_accuracy_af_af_resource.tsv"), Arrays.asList("--resource", TEST_DIR + "af_resource.vcf"), Collections.emptyMap()},
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth_new_sample_names.vcf", new File(TEST_DIR + "expected_correlations_basic.tsv"), new File(TEST_DIR + "expected_accuracy_basic.tsv"),
                        expectedAccuracyAFBasic, Arrays.asList("--sample-map", "SAMPLE1:SAMPLE3", "--sample-map", "SAMPLE2:SAMPLE4"), Collections.emptyMap()},
                {TEST_DIR + "eval_single_sample.vcf", TEST_DIR + "truth_single_sample.vcf", new File(TEST_DIR + "expected_correlations_single_sample.tsv"), new File(TEST_DIR + "expected_accuracy_single_sample.tsv"),
                        new File(TEST_DIR + "expected_accuracy_af_single_sample.tsv"), Arrays.asList("--force-compare-single-sample"), Collections.emptyMap()},
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth_multiallelic.vcf", new File(TEST_DIR + "expected_correlations_basic.tsv"), new File(TEST_DIR + "expected_accuracy_basic.tsv"),
                        expectedAccuracyAFBasic, Collections.emptyList(), Collections.emptyMap()}, //multiallelic truth should give the same results as split truth
                {TEST_DIR + "eval_monomorphic.vcf", TEST_DIR + "truth.vcf", new File(TEST_DIR + "expected_correlations_af_resource.tsv"), new File(TEST_DIR + "expected_accuracy_basic.tsv"),
                        new File(TEST_DIR + "expected_accuracy_af_af_resource.tsv"), Arrays.asList("--resource", TEST_DIR + "af_resource.vcf"), Collections.emptyMap()}, //eval records with no alt allele should be skipped
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth.vcf", new File(TEST_DIR + "expected_correlations_basic.tsv"), new File(TEST_DIR + "expected_accuracy_basic.tsv"),
                        expectedAccuracyAFBasic, Arrays.asList("--ac-resource", TEST_DIR + "ac_resource.vcf", "--max-ac-bin", "3"),
                        Collections.singletonMap("output-accuracy-ac", new File(TEST_DIR + "expected_accuracy_ac.tsv"))},
                {TEST_DIR + "eval_gp.vcf", TEST_DIR + "truth.vcf", new File(TEST_DIR + "expected_correlations_basic.tsv"), new File(TEST_DIR + "expected_accuracy_basic.tsv"),
                        expectedAccuracyAFBasic, Arrays.asList("--n-calibration-bins", "2"),
                        Collections.singletonMap("output-gp-calibration", new File(TEST_DIR + "expected_gp_calibration.tsv"))},
                {TEST_DIR + "eval_ids.vcf", TEST_DIR + "truth.vcf", new File(TEST_DIR + "expected_correlations_keep_ids.tsv"), new File(TEST_DIR + "expected_accuracy_keep_ids.tsv"),
                        new File(TEST_DIR + "expected_accuracy_af_keep_ids.tsv"),
                        Arrays.asList("--keep-ids", "rs100", "--keep-ids", "rs105", "--keep-ids", "rs415", "--keep-ids", "rs500"), Collections.emptyMap()}
        };
    }

    /**
     * @param additionalOutputs map from the long name of an optional output argument to the corresponding expected output
     */
    @Test(dataProvider = "integrationTestDataProvider")
    void integrationTests(final String eval, final String truth, final File expectedCorrelation, final File expectedAccuracy, final File expectedAccuracyAF,
                          final List<String> additionalArgs, final Map<String, File> additionalOutputs) throws IOException {
        final File outputCorrelation = createTempFile("correlation", ".tsv");
        final File outputAccuracy = createTempFile("accuracy", ".tsv");
        final File outputAccuracyAF = createTempFile("accuracy_af", ".tsv");

        final List<String> args = new ArrayList<>(additionalArgs);
        final Map<File, File> additionalOutputFiles = new LinkedHashMap<>();
        for (final Map.Entry<String, File> entry : additionalOutputs.entrySet()) {
            final File output = createTempFile(entry.getKey(), ".tsv");
            args.addAll(Arrays.asList("--" + entry.getKey(), output.getAbsolutePath()));
            additionalOutputFiles.put(output, entry.getValue());
        }

        runCommandLine(buildArgs(eval, truth, outputCorrelation.getAbsolutePath(), outputAccuracy.getAbsolutePath(), outputAccuracyAF.getAbsolutePath(), args));

        IntegrationTestSpec.assertEqualTextFiles(outputCorrelation, expectedCorrelation, "#");
        IntegrationTestSpec.assertEqualTextFiles(outputAccuracy, expectedAccuracy, "#");
        IntegrationTestSpec.assertEqualTextFiles(outputAccuracyAF, expectedAccuracyAF, "#");
        for (final Map.Entry<File, File> entry : additionalOutputFiles.entrySet()) {
            IntegrationTestSpec.assertEqualTextFiles(entry.getKey(), entry.getValue(), "#");
        }
    }


    @DataProvider(name = "exceptionProducingDataProvider")
    Object[][] exceptionProductionDataProvider() {
        return new Object[][] {
                {TEST_DIR + "eval_X.vcf", TEST_DIR + "truth_X.vcf", Collections.emptyList()}, //differing ploidy without argument to allow
                {TEST_DIR + "eval_single_sample.vcf", TEST_DIR + "truth_single_sample.vcf", Collections.emptyList()}, //single sample differing samples without argument to allow
                {TEST_DIR + "eval_single_sample.vcf", TEST_DIR + "truth.vcf", Arrays.asList("--force-compare-single-sample")}, //single sample argument with multi-sample truth
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth_single_sample.vcf", Arrays.asList("--force-compare-single-sample")}, //single sample argument with multi-sample eval
                {TEST_DIR + "eval_single_sample.vcf", TEST_DIR + "truth_single_sample.vcf", Arrays.asList("--force-compare-single-sample", "--sample-map", "SAMPLE1:SAMPLE3")}, //single sample argument with sample mapping
                {TEST_DIR + "eval.vcf", TEST_DIR + "truth_single_sample.vcf", Collections.emptyList()} //no comparisons to be run
        };
    }
    @Test(dataProvider = "exceptionProducingDataProvider", expectedExceptions = GATKException.class)
    void exceptionProducingTests(final String eval, final String truth, final List<String> additionalArgs) {
        runCommandLine(buildArgs(eval, truth,
                createTempFile("correlation", ".tsv").getAbsolutePath(), createTempFile("accuracy", ".tsv").getAbsolutePath(), createTempFile("accuracy_af", ".tsv").getAbsolutePath(),
                additionalArgs));
    }

    List<String> buildArgs(final String evalPath, final String truthPath, final String outputCorrelationPath, final String outputAccuracyPath,
                           final String outputAccuracyAFPath, final List<String> additionalArgs) {
        final List<String> args = new ArrayList<>(Arrays.asList(
                "--first-bin-right-edge", "0.1",
                "-nbins", "5",
                "--af-annotations", "SAMPLE1:AF1",
                "--af-annotations", "SAMPLE2:AF2",
                "--eval", evalPath,
                "--truth", truthPath,
                "--O", outputCorrelationPath,
                "--OA", outputAccuracyPath,
                "--output-accuracy-af", outputAccuracyAFPath
        ));
        args.addAll(additionalArgs);

        return args;
    }
}
