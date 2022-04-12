package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.collect.Lists;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.VariantType;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collectors;

/**
 * See documentation for {@link ExtractVariantAnnotationsIntegrationTest} for information about how inputs and
 * expected outputs used there are related to those used here and in {@link TrainVariantAnnotationsModelIntegrationTest}.
 */
public final class ScoreVariantAnnotationsIntegrationTest extends CommandLineProgramTest {

    // If true, update the expected outputs in tests that assert an exact match vs. prior output,
    // instead of actually running the tests. Can be used with "./gradlew test -Dtest.single=ScoreVariantAnnotationsIntegrationTest"
    // to update all of the exact-match tests at once. After you do this, you should look at the
    // diffs in the new expected outputs in git to confirm that they are consistent with expectations.
    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    /**
     * Make sure that someone didn't leave the UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS toggle turned on.
     */
    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled.");
    }

    private static final File PACKAGE_TEST_FILES_DIR = new File(largeFileTestDir,
            "org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/");
    private static final File TEST_FILES_DIR = new File(largeFileTestDir,
            "org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/score");
    private static final File INPUT_FROM_TRAIN_EXPECTED_TEST_FILES_DIR = new File(largeFileTestDir,
            "org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/train/expected");
    private static final File EXPECTED_TEST_FILES_DIR = new File(TEST_FILES_DIR, "expected");

    private static final File ISOLATION_FOREST_PYTHON_SCRIPT = new File(packageMainResourcesDir,
            "tools/walkers/vqsr/scalable/isolation-forest.py");

    private static final File INPUT_VCF = new File(PACKAGE_TEST_FILES_DIR, "input/small_callset_low_threshold.sites-only.chr1.1-10M.vcf");

    // Supplier and functions for creating and adding various arguments to an ArgumentsBuilder.
    private static final Supplier<ArgumentsBuilder> BASE_ARGUMENTS_BUILDER_SUPPLIER = () -> {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addVCF(INPUT_VCF);
        argsBuilder.addFlag(LabeledVariantAnnotationsWalker.DO_NOT_GZIP_VCF_OUTPUT_LONG_NAME);
        argsBuilder.add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);
        return argsBuilder;
    };
    private static final BiFunction<ArgumentsBuilder, String, ArgumentsBuilder> ADD_MODEL_PREFIX = (argsBuilder, modelPrefix) -> {
        argsBuilder.add(ScoreVariantAnnotations.MODEL_PREFIX_LONG_NAME, modelPrefix);
        return argsBuilder;
    };
    private static final BiFunction<ArgumentsBuilder, Double, ArgumentsBuilder> ADD_CALIBRATION_SENSITIVITY_THRESHOLD = (argsBuilder, calibrationSensitivityThreshold) -> {
        argsBuilder.add(ScoreVariantAnnotations.CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME, calibrationSensitivityThreshold);
        return argsBuilder;
    };
    private static final Function<ArgumentsBuilder, ArgumentsBuilder> ADD_ISOLATION_FOREST_PYTHON_SCRIPT = argsBuilder -> {
        argsBuilder.add(ScoreVariantAnnotations.PYTHON_SCRIPT_LONG_NAME, ISOLATION_FOREST_PYTHON_SCRIPT);
        return argsBuilder;
    };

    /**
     * Exact-match tests for (non-exhaustive) configurations given by the Cartesian product of the following options:
     * 1) non-allele-specific vs. allele-specific
     * 2) Java Bayesian Gaussian Mixture Model (BGMM) backend vs. python sklearn IsolationForest backend
     * 3) SNP-only vs. SNP+INDEL (for both of these options, we use trained models that contain both SNP and INDEL scorers as input)
     *  TODO the BGMM has been reduced to a stub for this initial PR; subsequent PRs will cover the backend code and reconnect the stub
     *  TODO warm-start BGMM?
     */
    @DataProvider(name = "dataValidInputs")
    public Object[][] dataValidInputs() {
        final List<List<Pair<String, Function<ArgumentsBuilder, ArgumentsBuilder>>>> testConfigurations = Lists.cartesianProduct(
                Arrays.asList(
                        Pair.of("extract.nonAS.snpIndel.posUn.train.snpIndel.posNeg.IF.score", ADD_ISOLATION_FOREST_PYTHON_SCRIPT), // TODO add BGMM
                        Pair.of("extract.AS.snpIndel.posUn.train.snpIndel.posNeg.IF.score", ADD_ISOLATION_FOREST_PYTHON_SCRIPT)),
                Arrays.asList(
                        Pair.of("snp", ExtractVariantAnnotationsIntegrationTest.ADD_SNP_MODE_AND_RESOURCES),
                        Pair.of("snpIndel", ExtractVariantAnnotationsIntegrationTest.ADD_SNP_MODE_AND_RESOURCES
                                .andThen(ExtractVariantAnnotationsIntegrationTest.ADD_INDEL_MODE_AND_RESOURCES))));

        return testConfigurations.stream()
                .map(tagAndAddFunctionPairs -> new Object[]{
                        tagAndAddFunctionPairs.stream().map(Pair::getLeft).collect(Collectors.joining(".")), // e.g., extract.nonAS.snpIndel.posUn.train.snpIndel.posNeg.IF.score.snp
                        tagAndAddFunctionPairs.stream().map(Pair::getRight)                                              // creates the corresponding ArgumentsBuilder
                                .reduce(Function.identity(), Function::andThen)                                          //  by stringing together functions that add the
                                .apply(BASE_ARGUMENTS_BUILDER_SUPPLIER.get())})                                          //  appropriate arguments
                .toArray(Object[][]::new);
    }

    @Test(dataProvider = "dataValidInputs", groups = {"python"}) // python environment is required to run tool and to use h5diff for exact-match comparisons
    public void testValidInputs(final String tag,
                                final ArgumentsBuilder argsBuilder) {
        final File outputDir = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? EXPECTED_TEST_FILES_DIR : createTempDir("score");
        final String outputPrefix = String.format("%s/%s", outputDir, tag);
        argsBuilder.addOutput(outputPrefix);

        // add arguments for model prefix based on the
        // train tag (the portion of the tag preceding ".score", e.g., extract.nonAS.snpIndel.posUn.train.posNeg.IF),
        // which gives the basename for the model files
        final String trainTag = tag.split(".score")[0];
        if (tag.contains("nonAS")) {
            ExtractVariantAnnotationsIntegrationTest.ADD_NON_ALLELE_SPECIFIC_ANNOTATIONS.apply(argsBuilder);
        } else {
            ExtractVariantAnnotationsIntegrationTest.ADD_ALLELE_SPECIFIC_ANNOTATIONS.apply(argsBuilder);
        }
        final String modelPrefix = new File(INPUT_FROM_TRAIN_EXPECTED_TEST_FILES_DIR, trainTag).toString();
        final Function<ArgumentsBuilder, ArgumentsBuilder> addModelPrefix = ab ->
                ADD_MODEL_PREFIX.apply(ab, modelPrefix);
        final double calibrationSensitivityThreshold = 0.9;
        final Function<ArgumentsBuilder, ArgumentsBuilder> addCalibrationSensitivityThreshold = ab ->
                ADD_CALIBRATION_SENSITIVITY_THRESHOLD.apply(ab, calibrationSensitivityThreshold);
        addModelPrefix.andThen(addCalibrationSensitivityThreshold).apply(argsBuilder);

        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, "INFO");
        runCommandLine(argsBuilder);

        if (!UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            assertExpectedOutputs(tag, outputPrefix);
        }
    }

    private static void assertExpectedOutputs(final String tag,
                                              final String outputPrefix) {
        SystemCommandUtilsTest.runSystemCommand(String.format("diff %s/%s.vcf %s.vcf",
                EXPECTED_TEST_FILES_DIR, tag, outputPrefix));
//        SystemCommandUtilsTest.runSystemCommand(String.format("diff %s/%s.vcf.idx %s.vcf.idx",
//                EXPECTED_TEST_FILES_DIR, tag, outputPrefix));
        SystemCommandUtilsTest.runSystemCommand(String.format("h5diff %s/%s.annot.hdf5 %s.annot.hdf5",
                EXPECTED_TEST_FILES_DIR, tag, outputPrefix));
        SystemCommandUtilsTest.runSystemCommand(String.format("h5diff %s/%s.scores.hdf5 %s.scores.hdf5",
                EXPECTED_TEST_FILES_DIR, tag, outputPrefix));
    }
}