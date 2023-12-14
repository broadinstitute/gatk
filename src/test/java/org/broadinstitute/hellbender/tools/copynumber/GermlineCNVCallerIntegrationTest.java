package org.broadinstitute.hellbender.tools.copynumber;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.CopyNumberTestUtils;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.arguments.GermlineDenoisingModelArgumentCollection;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;


/**
 * Integration tests for {@link GermlineCNVCaller}.
 */
public final class GermlineCNVCallerIntegrationTest extends CommandLineProgramTest {

    // If true, update the expected outputs in tests that assert an exact or approximate match vs. prior output,
    // instead of actually running the tests.
    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    private static final String GCNV_SIM_DATA_DIR = toolsTestDir + "copynumber/gcnv-sim-data/";
    private static final File[] TEST_COUNT_FILES = IntStream.range(0, 20)
            .mapToObj(n -> new File(GCNV_SIM_DATA_DIR + String.format("SAMPLE_%03d_counts.tsv", n)))
            .toArray(File[]::new);
    private static final File CONTIG_PLOIDY_CALLS_OUTPUT_DIR = new File(GCNV_SIM_DATA_DIR + "contig-ploidy-calls/");
    private static final File SIM_INTERVAL_LIST_SUBSET_FILE = new File(GCNV_SIM_DATA_DIR + "sim_intervals_subset.interval_list");
    private static final File OUTPUT_DIR = createTempDir("test-germline-cnv");
    private static final String SAMPLE_DIR_EXT = "SAMPLE_0";
    private static final String GCNV_TEST_OUTPUT_DIR = toolsTestDir + "copynumber/gcnv-numerical-accuracy/";
    private static final String TEST_SHARD = "shard_0";
    private static final File MODEL_EXACT_MATCH_EXPECTED_SUB_DIR = new File(GCNV_TEST_OUTPUT_DIR + TEST_SHARD + "-model/");
    private static final File CALLS_EXACT_MATCH_EXPECTED_SUB_DIR = new File(GCNV_TEST_OUTPUT_DIR + TEST_SHARD + "-calls/" + SAMPLE_DIR_EXT);
    private static final File SIM_INTERVAL_LIST_SHARD_0 = new File(GCNV_SIM_DATA_DIR + "sim_intervals_shard_0.interval_list");
    private static final File SIM_INTERVAL_LIST_SHARD_0_ANNOTATED_FILE = new File(GCNV_SIM_DATA_DIR + "sim_intervals_shard_0.annotated.tsv");
    private static final double ALLOWED_DELTA_FOR_DOUBLE_VALUES = 1E-6;

    final List<String> MODEL_FILES_TO_COMPARE = Arrays.asList("log_q_tau_tk.tsv", "mu_ard_u_interval__.tsv", "mu_psi_t_log__.tsv",
            "std_ard_u_interval__.tsv", "std_psi_t_log__.tsv", "mu_W_tu.tsv", "mu_log_mean_bias_t.tsv", "std_W_tu.tsv", "std_log_mean_bias_t.tsv");
    final List<String> CALLS_FILES_TO_COMPARE = Arrays.asList("baseline_copy_number_t.tsv", "log_c_emission_tc.tsv",
            "log_q_c_tc.tsv", "mu_psi_s_log__.tsv", "mu_read_depth_s_log__.tsv",
            "mu_z_su.tsv", "sample_name.txt", "std_psi_s_log__.tsv",
            "std_read_depth_s_log__.tsv", "std_z_su.tsv", "mu_denoised_copy_ratio_t.tsv", "std_denoised_copy_ratio_t.tsv");

    /**
     * Make sure that someone didn't leave the {@value UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS} toggle turned on
     */
    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    /**
     * Run the tool in the COHORT mode for all 20 samples on a small subset of intervals
     */
    @Test(groups = {"python"})
    public void testCohortWithoutIntervalAnnotations() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        Arrays.stream(TEST_COUNT_FILES).forEach(argsBuilder::addInput);
        argsBuilder.add(GermlineCNVCaller.RUN_MODE_LONG_NAME, GermlineCNVCaller.RunMode.COHORT.name())
                .add("L", SIM_INTERVAL_LIST_SUBSET_FILE)
                .add(GermlineCNVCaller.CONTIG_PLOIDY_CALLS_DIRECTORY_LONG_NAME,
                        CONTIG_PLOIDY_CALLS_OUTPUT_DIR.getAbsolutePath())
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, OUTPUT_DIR.getAbsolutePath())
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test-germline-cnv-cohort")
                .add(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString());
        runCommandLine(argsBuilder);
    }

    /**
     * Run the tool in the COHORT mode for a case when number of provided intervals is two (minimum number allowed).
     */
    @Test(groups = {"python"})
    public void testCohortWithTwoIntervals() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        Arrays.stream(TEST_COUNT_FILES, 0, 5).forEach(argsBuilder::addInput);
        argsBuilder.add(GermlineCNVCaller.RUN_MODE_LONG_NAME, GermlineCNVCaller.RunMode.COHORT.name())
                .add("L", "1:620998-622131")
                .add("L", "1:877712-877902")
                .add(GermlineCNVCaller.CONTIG_PLOIDY_CALLS_DIRECTORY_LONG_NAME,
                        CONTIG_PLOIDY_CALLS_OUTPUT_DIR.getAbsolutePath())
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, OUTPUT_DIR.getAbsolutePath())
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test-germline-cnv-cohort-two-targets")
                .add(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString());
        runCommandLine(argsBuilder);
    }

    /**
     * Run the tool in CASE mode for the first 5 samples using the model generated by
     * {@link #testCohortWithoutIntervalAnnotations()}
     */
    @Test(groups = {"python"}, dependsOnMethods = "testCohortWithoutIntervalAnnotations")
    public void testCase() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        Arrays.stream(TEST_COUNT_FILES, 0, 5).forEach(argsBuilder::addInput);
        argsBuilder.add(GermlineCNVCaller.RUN_MODE_LONG_NAME, GermlineCNVCaller.RunMode.CASE.name())
                .add(GermlineCNVCaller.CONTIG_PLOIDY_CALLS_DIRECTORY_LONG_NAME,
                        CONTIG_PLOIDY_CALLS_OUTPUT_DIR.getAbsolutePath())
                .add(CopyNumberStandardArgument.MODEL_LONG_NAME,
                        new File(OUTPUT_DIR, "test-germline-cnv-cohort-model").getAbsolutePath())
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, OUTPUT_DIR.getAbsolutePath())
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test-germline-cnv-case");
        runCommandLine(argsBuilder);
    }

    @Test(groups = {"python"}, expectedExceptions = IllegalArgumentException.class)
    public void testCaseWithoutModel() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        Arrays.stream(TEST_COUNT_FILES, 0, 5).forEach(argsBuilder::addInput);
        argsBuilder.add(GermlineCNVCaller.RUN_MODE_LONG_NAME, GermlineCNVCaller.RunMode.CASE.name())
                .add(GermlineCNVCaller.CONTIG_PLOIDY_CALLS_DIRECTORY_LONG_NAME,
                        CONTIG_PLOIDY_CALLS_OUTPUT_DIR.getAbsolutePath())
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, OUTPUT_DIR.getAbsolutePath())
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test-germline-cnv-case");
        runCommandLine(argsBuilder);
    }

    @Test(groups = {"python"}, expectedExceptions = IllegalArgumentException.class)
    public void testCaseWithHiddenArguments() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        Arrays.stream(TEST_COUNT_FILES, 0, 5).forEach(argsBuilder::addInput);
        argsBuilder.add(GermlineCNVCaller.RUN_MODE_LONG_NAME, GermlineCNVCaller.RunMode.CASE.name())
                .add(GermlineCNVCaller.CONTIG_PLOIDY_CALLS_DIRECTORY_LONG_NAME,
                        CONTIG_PLOIDY_CALLS_OUTPUT_DIR.getAbsolutePath())
                .add(CopyNumberStandardArgument.MODEL_LONG_NAME,
                        new File(OUTPUT_DIR, "test-germline-cnv-cohort-model").getAbsolutePath())
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, OUTPUT_DIR.getAbsolutePath())
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test-germline-cnv-case");
        // add argument that is not applicable in CASE mode
        argsBuilder.add(GermlineDenoisingModelArgumentCollection.INTERVAL_PSI_SCALE_LONG_NAME, 0.1);
        runCommandLine(argsBuilder);
    }

    /**
     * Test that {@link GermlineCNVCaller} outputs fall within a numerical accuracy bound defined. This test is meant
     * to detect things like major Python library updates that affect gCNV results, and detect unintentional
     * consequences of minor model changes.
     */
    @Test(groups = {"python"})
    public void testNumericalAccuracy() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        final String outputPrefix = TEST_SHARD;
        Arrays.stream(TEST_COUNT_FILES).forEach(argsBuilder::addInput);
        final String outputDir = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? GCNV_TEST_OUTPUT_DIR : OUTPUT_DIR.getAbsolutePath();

        argsBuilder.add(GermlineCNVCaller.RUN_MODE_LONG_NAME, GermlineCNVCaller.RunMode.COHORT.name())
                .add(StandardArgumentDefinitions.INTERVALS_SHORT_NAME, SIM_INTERVAL_LIST_SHARD_0)
                .add(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString())
                .add(GermlineCNVCaller.CONTIG_PLOIDY_CALLS_DIRECTORY_LONG_NAME,
                        CONTIG_PLOIDY_CALLS_OUTPUT_DIR.getAbsolutePath())
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, outputDir)
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, outputPrefix);
        runCommandLine(argsBuilder);

        // Test that values of outputs are approximately numerically equivalent
        if (!UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            MODEL_FILES_TO_COMPARE.forEach(f -> CopyNumberTestUtils.assertFilesEqualUpToAllowedDeltaForDoubleValues(
                    new File(Paths.get(OUTPUT_DIR.getAbsolutePath(), outputPrefix + "-model").toString(), f),
                    new File(MODEL_EXACT_MATCH_EXPECTED_SUB_DIR, f),
                    ALLOWED_DELTA_FOR_DOUBLE_VALUES,
                    logger));
            CALLS_FILES_TO_COMPARE.forEach(f -> CopyNumberTestUtils.assertFilesEqualUpToAllowedDeltaForDoubleValues(
                    new File(Paths.get(OUTPUT_DIR.getAbsolutePath(), outputPrefix + "-calls", SAMPLE_DIR_EXT).toString(), f),
                    new File(CALLS_EXACT_MATCH_EXPECTED_SUB_DIR, f),
                    ALLOWED_DELTA_FOR_DOUBLE_VALUES,
                    logger));
        } else {  
            // remove files not used for testing
            try {
                FileUtils.deleteDirectory(new File(Paths.get(GCNV_TEST_OUTPUT_DIR, outputPrefix + "-tracking").toString()));
            } catch (final IOException ex) {
                throw new GATKException("Could not remove GermlineCNVCaller tracking files.");
            }
            IntStream.range(1, 20).forEach(
                    s -> {
                        try {
                            FileUtils.deleteDirectory(new File(Paths.get(GCNV_TEST_OUTPUT_DIR, outputPrefix + "-calls", "SAMPLE_" + s).toString()));
                        } catch (final IOException ex) {
                            throw new GATKException(String.format("Could not remove a GermlineCNVCaller call directory for sample %d", s));
                        }
                    });
        }
    }

    /**
     * Run the tool in the COHORT mode using a provided model for initialization. Note for developers: the model used
     * for initialization is created by the {@link #testCohortWithoutIntervalAnnotations()} test.
     */
    @Test(groups = {"python"}, dependsOnMethods = "testCohortWithoutIntervalAnnotations")
    public void testCohortWithInputModel() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        Arrays.stream(TEST_COUNT_FILES).forEach(argsBuilder::addInput);
        argsBuilder.add(GermlineCNVCaller.RUN_MODE_LONG_NAME, GermlineCNVCaller.RunMode.COHORT.name())
                .add(StandardArgumentDefinitions.INTERVALS_SHORT_NAME, SIM_INTERVAL_LIST_SUBSET_FILE)
                .add(CopyNumberStandardArgument.MODEL_LONG_NAME,
                        new File(OUTPUT_DIR, "test-germline-cnv-cohort-model").getAbsolutePath())
                .add(GermlineCNVCaller.CONTIG_PLOIDY_CALLS_DIRECTORY_LONG_NAME,
                        CONTIG_PLOIDY_CALLS_OUTPUT_DIR.getAbsolutePath())
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, OUTPUT_DIR.getAbsolutePath())
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test-germline-cnv-cohort-with-input-model")
                .add(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString());
        runCommandLine(argsBuilder);
    }

    @Test(groups = {"python"})
    public void testCohortWithAnnotatedIntervals() {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        Arrays.stream(TEST_COUNT_FILES).forEach(argsBuilder::addInput);
        argsBuilder.add(GermlineCNVCaller.RUN_MODE_LONG_NAME, GermlineCNVCaller.RunMode.COHORT.name())
                .add(StandardArgumentDefinitions.INTERVALS_SHORT_NAME, SIM_INTERVAL_LIST_SHARD_0)
                .add(CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_LONG_NAME, SIM_INTERVAL_LIST_SHARD_0_ANNOTATED_FILE)
                .add(GermlineCNVCaller.CONTIG_PLOIDY_CALLS_DIRECTORY_LONG_NAME,
                        CONTIG_PLOIDY_CALLS_OUTPUT_DIR.getAbsolutePath())
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, OUTPUT_DIR.getAbsolutePath())
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test-germline-cnv-cohort")
                .add(IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString());
        runCommandLine(argsBuilder);
    }
}
