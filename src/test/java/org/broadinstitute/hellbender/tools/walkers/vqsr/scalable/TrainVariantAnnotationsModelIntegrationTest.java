package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.collect.Lists;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
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
 * expected outputs used there are related to those used here and in {@link ScoreVariantAnnotationsIntegrationTest}.
 */
public final class TrainVariantAnnotationsModelIntegrationTest extends CommandLineProgramTest {

    // If true, update the expected outputs in tests that assert an exact match vs. prior output,
    // instead of actually running the tests. Can be used with "./gradlew test -Dtest.single=TrainVariantAnnotationsIntegrationTest"
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

    private static final File TEST_FILES_DIR = new File(largeFileTestDir,
            "org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/train");
    private static final File INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR = new File(largeFileTestDir,
            "org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/extract/expected");
    private static final File EXPECTED_TEST_FILES_DIR = new File(TEST_FILES_DIR, "expected");

    private static final File ISOLATION_FOREST_PYTHON_SCRIPT = new File(packageMainResourcesDir,
            "tools/walkers/vqsr/scalable/isolation-forest.py");
    private static final File ISOLATION_FOREST_HYPERPARAMETERS_JSON = new File(TEST_FILES_DIR,
            "isolation-forest-hyperparameters.json");

    // Supplier and functions for creating and adding various arguments to an ArgumentsBuilder.
    private static final Supplier<ArgumentsBuilder> BASE_ARGUMENTS_BUILDER_SUPPLIER = ArgumentsBuilder::new;

    private static final BiFunction<ArgumentsBuilder, File, ArgumentsBuilder> ADD_ANNOTATIONS_HDF5 = (argsBuilder, annotationsHDF5) -> {
        argsBuilder.add(TrainVariantAnnotationsModel.ANNOTATIONS_HDF5_LONG_NAME, annotationsHDF5);
        return argsBuilder;
    };
    private static final BiFunction<ArgumentsBuilder, File, ArgumentsBuilder> ADD_UNLABELED_ANNOTATIONS_HDF5 = (argsBuilder, unlabeledAnnotationsHDF5) -> {
        argsBuilder.add(TrainVariantAnnotationsModel.UNLABELED_ANNOTATIONS_HDF5_LONG_NAME, unlabeledAnnotationsHDF5);
        return argsBuilder;
    };
    private static final BiFunction<ArgumentsBuilder, Double, ArgumentsBuilder> ADD_CALIBRATION_SENSITIVITY_THRESHOLD = (argsBuilder, calibrationSensitivityThreshold) -> {
        argsBuilder.add(TrainVariantAnnotationsModel.CALIBRATION_SENSITIVITY_THRESHOLD_LONG_NAME, calibrationSensitivityThreshold);
        return argsBuilder;
    };
    private static final Function<ArgumentsBuilder, ArgumentsBuilder> ADD_SNP_MODE = argsBuilder -> {
        argsBuilder.add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, VariantType.SNP);
        return argsBuilder;
    };
    private static final Function<ArgumentsBuilder, ArgumentsBuilder> ADD_INDEL_MODE = argsBuilder -> {
        argsBuilder.add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, VariantType.INDEL);
        return argsBuilder;
    };
    private static final Function<ArgumentsBuilder, ArgumentsBuilder> ADD_ISOLATION_FOREST_PYTHON_SCRIPT = argsBuilder -> {
        argsBuilder.add(TrainVariantAnnotationsModel.PYTHON_SCRIPT_LONG_NAME, ISOLATION_FOREST_PYTHON_SCRIPT);
        return argsBuilder;
    };
    private static final Function<ArgumentsBuilder, ArgumentsBuilder> ADD_ISOLATION_FOREST_HYPERPARAMETERS_JSON = argsBuilder -> {
        argsBuilder.add(TrainVariantAnnotationsModel.HYPERPARAMETERS_JSON_LONG_NAME, ISOLATION_FOREST_HYPERPARAMETERS_JSON);
        return argsBuilder;
    };

    /**
     * Exact-match tests for (non-exhaustive) configurations given by the Cartesian product of the following options:
     *  1) non-allele-specific vs. allele-specific
     *  2) SNP-only vs. SNP+INDEL (for both of these options, we use extracted annotations that contain both SNP and INDEL variants as input)
     *  3) positive (training with {extract-tag}.annot.hdf5) vs. positive-unlabeled (training with {extract-tag}.annot.hdf5 and {extract-tag}.unlabeled.annot.hdf5)
     *  4) Java Bayesian Gaussian Mixture Model (BGMM) backend vs. python sklearn IsolationForest backend
     *  TODO the BGMM has been reduced to a stub for this initial PR; subsequent PRs will cover the backend code and reconnect the stub
     *  TODO warm-start BGMM?
     */
    @DataProvider(name = "dataValidInputs")
    public Object[][] dataValidInputs() {
        final List<List<Pair<String, Function<ArgumentsBuilder, ArgumentsBuilder>>>> testConfigurations = Lists.cartesianProduct(
                Arrays.asList(
                        Pair.of("extract.nonAS.snpIndel.posUn.train", Function.identity()),
                        Pair.of("extract.nonAS.snpIndel.posUn.train", Function.identity()),
                        Pair.of("extract.AS.snpIndel.posUn.train", Function.identity()),
                        Pair.of("extract.AS.snpIndel.posUn.train", Function.identity())),
                Arrays.asList(
                        Pair.of("snp", ADD_SNP_MODE),
                        Pair.of("snpIndel", ADD_SNP_MODE.andThen(ADD_INDEL_MODE))),
                Arrays.asList(              // we will consume the tag and add appropriate arguments for positive and positive-negative training below
                        Pair.of("posOnly", Function.identity()),
                        Pair.of("posNeg", Function.identity())),
                Collections.singletonList(
                        Pair.of("IF", ADD_ISOLATION_FOREST_PYTHON_SCRIPT.andThen(ADD_ISOLATION_FOREST_HYPERPARAMETERS_JSON))));

        return testConfigurations.stream()
                .map(tagAndAddFunctionPairs -> new Object[]{
                        tagAndAddFunctionPairs.stream().map(Pair::getLeft).collect(Collectors.joining(".")), // e.g., extract.nonAS.snpIndel.posUn.train.snp.posOnly.IF
                        tagAndAddFunctionPairs.stream().map(Pair::getRight)                                              // creates the corresponding ArgumentsBuilder
                                .reduce(Function.identity(), Function::andThen)                                          //  by stringing together functions that add the
                                .apply(BASE_ARGUMENTS_BUILDER_SUPPLIER.get())})                                               //  appropriate arguments
                .toArray(Object[][]::new);
    }

    @Test(dataProvider = "dataValidInputs")
    public void testValidInputs(final String tag,
                                final ArgumentsBuilder argsBuilder) {
        final File outputDir = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? EXPECTED_TEST_FILES_DIR : createTempDir("train");
        final String outputPrefix = String.format("%s/%s", outputDir, tag);
        argsBuilder.addOutput(outputPrefix);

        // add arguments for positive/unlabeled annotations based on the
        // extract tag (the portion of the tag preceding ".train", e.g., extract.nonAS.snpIndel.posUn),
        // which gives the basename for the annotation files
        final String extractTag = tag.split(".train")[0];
        final File positiveAnnotationsHDF5 = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                String.format("%s.annot.hdf5", extractTag));
        final Function<ArgumentsBuilder, ArgumentsBuilder> addPositiveAnnotations = ab ->
                ADD_ANNOTATIONS_HDF5.apply(ab, positiveAnnotationsHDF5);
        if (tag.contains("posNeg")) {
            final File unlabeledAnnotationsHDF5 = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                    String.format("%s.unlabeled.annot.hdf5", extractTag));
            final Function<ArgumentsBuilder, ArgumentsBuilder> addUnlabeledAnnotations = ab ->
                    ADD_UNLABELED_ANNOTATIONS_HDF5.apply(ab, unlabeledAnnotationsHDF5);
            final double calibrationSensitivityThreshold = 0.9;
            final Function<ArgumentsBuilder, ArgumentsBuilder> addCalibrationSensitivityThreshold = ab ->
                    ADD_CALIBRATION_SENSITIVITY_THRESHOLD.apply(ab, calibrationSensitivityThreshold);
            addPositiveAnnotations.andThen(addUnlabeledAnnotations).andThen(addCalibrationSensitivityThreshold).apply(argsBuilder);
        } else {
            addPositiveAnnotations.apply(argsBuilder);
        }

        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, "INFO");
        runCommandLine(argsBuilder);

        if (!UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            assertExpectedOutputs(tag, outputPrefix);
        }
    }

    private static void assertExpectedOutputs(final String tag,
                                              final String outputPrefix) {
        if (tag.contains("train.snp")) {
            assertExpectedOutputsForVariantType(tag, outputPrefix, "snp");
            assertOutputsForVariantTypeDoNotExist(outputPrefix, "indel");
        } else if (tag.contains("train.snpIndel")) {
            assertExpectedOutputsForVariantType(tag, outputPrefix, "snp");
            assertExpectedOutputsForVariantType(tag, outputPrefix, "indel");
        } else {
            Assert.fail("Unknown variant-type tag.");
        }
    }

    private static void assertExpectedOutputsForVariantType(final String tag,
                                                            final String outputPrefix,
                                                            final String variantType) {
        final String tagAndVariantType = String.format("%s.%s", tag, variantType);
        final String outputPrefixAndVariantType = String.format("%s.%s", outputPrefix, variantType);

        SystemCommandUtilsTest.runSystemCommand(String.format("h5diff %s/%s.trainingScores.hdf5 %s.trainingScores.hdf5",
                EXPECTED_TEST_FILES_DIR, tagAndVariantType, outputPrefixAndVariantType));
        SystemCommandUtilsTest.runSystemCommand(String.format("h5diff %s/%s.calibrationScores.hdf5 %s.calibrationScores.hdf5",
                EXPECTED_TEST_FILES_DIR, tagAndVariantType, outputPrefixAndVariantType));

        assertScorerExpectedOutputs(tagAndVariantType, outputPrefixAndVariantType, false);

        if (tag.contains("posNeg")) {
            SystemCommandUtilsTest.runSystemCommand(String.format("h5diff %s/%s.unlabeledScores.hdf5 %s.unlabeledScores.hdf5",
                    EXPECTED_TEST_FILES_DIR, tagAndVariantType, outputPrefixAndVariantType));
            assertScorerExpectedOutputs(tagAndVariantType, outputPrefixAndVariantType, true);
        } else {
            Assert.assertFalse(new File(outputPrefixAndVariantType, ".unlabeledScores.hdf5").exists());
            Assert.assertFalse(new File(outputPrefixAndVariantType, ".negative.scorer.ser").exists());
            Assert.assertFalse(new File(outputPrefixAndVariantType, ".negative.scorer.pkl").exists());
        }
    }

    private static void assertOutputsForVariantTypeDoNotExist(final String outputPrefix,
                                                              final String variantType) {
        final String outputPrefixAndVariantType = String.format("%s.%s", outputPrefix, variantType);
        Assert.assertFalse(new File(outputPrefixAndVariantType, ".trainingScores.hdf5").exists());
        Assert.assertFalse(new File(outputPrefixAndVariantType, ".calibrationScores.hdf5").exists());
        Assert.assertFalse(new File(outputPrefixAndVariantType, ".unlabeledScores.hdf5").exists());
        Assert.assertFalse(new File(outputPrefixAndVariantType, ".scorer.ser").exists());
        Assert.assertFalse(new File(outputPrefixAndVariantType, ".scorer.pkl").exists());
        Assert.assertFalse(new File(outputPrefixAndVariantType, ".negative.scorer.ser").exists());
        Assert.assertFalse(new File(outputPrefixAndVariantType, ".negative.scorer.pkl").exists());
    }

    private static void assertScorerExpectedOutputs(final String tagAndVariantType,
                                                    final String outputPrefixAndVariantType,
                                                    final boolean isNegative) {
        final String positiveOrNegativeTag = isNegative ? ".negative" : "";
        if (tagAndVariantType.contains("BGMM")) {
            SystemCommandUtilsTest.runSystemCommand(String.format("diff %s/%s.scorer.ser %s.scorer.ser",
                    EXPECTED_TEST_FILES_DIR, tagAndVariantType + positiveOrNegativeTag, outputPrefixAndVariantType + positiveOrNegativeTag));
            Assert.assertFalse(new File(outputPrefixAndVariantType, ".scorer.pkl").exists());
        } else if (tagAndVariantType.contains("IF")) {
            SystemCommandUtilsTest.runSystemCommand(String.format("diff %s/%s.scorer.pkl %s.scorer.pkl",
                    EXPECTED_TEST_FILES_DIR, tagAndVariantType + positiveOrNegativeTag, outputPrefixAndVariantType + positiveOrNegativeTag));
            Assert.assertFalse(new File(outputPrefixAndVariantType, ".scorer.ser").exists());
        } else {
            Assert.fail("Unknown model-backend tag.");
        }
    }

    @Test
    public void testSNPOnlyModelsFromSNPOnlyAndSNPPlusIndelAnnotationsAreIdentical() {
        final File outputDir = createTempDir("train");

        final String outputPrefixSNPOnly = String.format("%s/test-snp", outputDir);
        final ArgumentsBuilder argsBuilderSNPOnly = BASE_ARGUMENTS_BUILDER_SUPPLIER.get();
        argsBuilderSNPOnly.addOutput(outputPrefixSNPOnly);
        final File positiveAnnotationsHDF5SNPOnly = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                "extract.nonAS.snp.pos.annot.hdf5");
        final Function<ArgumentsBuilder, ArgumentsBuilder> addPositiveAnnotationsSNPOnly = ab ->
                ADD_ANNOTATIONS_HDF5.apply(ab, positiveAnnotationsHDF5SNPOnly);
        ADD_ISOLATION_FOREST_PYTHON_SCRIPT
                .andThen(ADD_ISOLATION_FOREST_HYPERPARAMETERS_JSON)
                .andThen(addPositiveAnnotationsSNPOnly)
                .andThen(ADD_SNP_MODE)
                .apply(argsBuilderSNPOnly);
        runCommandLine(argsBuilderSNPOnly);

        final String outputPrefixSNPPlusIndel = String.format("%s/test-snpIndel", outputDir);
        final ArgumentsBuilder argsBuilderSNPPlusIndel = BASE_ARGUMENTS_BUILDER_SUPPLIER.get();
        argsBuilderSNPPlusIndel.addOutput(outputPrefixSNPPlusIndel);
        final File positiveAnnotationsHDF5SNPPlusIndel = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                "extract.nonAS.snpIndel.pos.annot.hdf5");
        final Function<ArgumentsBuilder, ArgumentsBuilder> addPositiveAnnotationsSNPPlusIndel = ab ->
                ADD_ANNOTATIONS_HDF5.apply(ab, positiveAnnotationsHDF5SNPPlusIndel);
        ADD_ISOLATION_FOREST_PYTHON_SCRIPT
                .andThen(ADD_ISOLATION_FOREST_HYPERPARAMETERS_JSON)
                .andThen(addPositiveAnnotationsSNPPlusIndel)
                .andThen(ADD_SNP_MODE)
                .apply(argsBuilderSNPPlusIndel);
        runCommandLine(argsBuilderSNPPlusIndel);

        SystemCommandUtilsTest.runSystemCommand(String.format("h5diff %s.snp.trainingScores.hdf5 %s.snp.trainingScores.hdf5",
                outputPrefixSNPOnly, outputPrefixSNPPlusIndel));
        SystemCommandUtilsTest.runSystemCommand(String.format("h5diff %s.snp.calibrationScores.hdf5 %s.snp.calibrationScores.hdf5",
                outputPrefixSNPOnly, outputPrefixSNPPlusIndel));
        SystemCommandUtilsTest.runSystemCommand(String.format("diff %s.snp.scorer.pkl %s.snp.scorer.pkl",
                outputPrefixSNPOnly, outputPrefixSNPPlusIndel));
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testUnlabeledAnnotationsSpecifiedWithoutCalibrationSensitivityThreshold() {
        final File outputDir = createTempDir("train");
        final String outputPrefix = String.format("%s/test", outputDir);
        final ArgumentsBuilder argsBuilder = BASE_ARGUMENTS_BUILDER_SUPPLIER.get();
        argsBuilder.addOutput(outputPrefix);
        final String extractTag = "extract.nonAS.snpIndel.posUn";
        final File positiveAnnotationsHDF5 = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                String.format("%s.annot.hdf5", extractTag));
        final Function<ArgumentsBuilder, ArgumentsBuilder> addPositiveAnnotations = ab ->
                ADD_ANNOTATIONS_HDF5.apply(ab, positiveAnnotationsHDF5);
        final File unlabeledAnnotationsHDF5 = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                String.format("%s.unlabeled.annot.hdf5", extractTag));
        final Function<ArgumentsBuilder, ArgumentsBuilder> addUnlabeledAnnotations = ab ->
                ADD_UNLABELED_ANNOTATIONS_HDF5.apply(ab, unlabeledAnnotationsHDF5);
        ADD_ISOLATION_FOREST_PYTHON_SCRIPT
                .andThen(ADD_ISOLATION_FOREST_HYPERPARAMETERS_JSON)
                .andThen(addPositiveAnnotations)
                .andThen(addUnlabeledAnnotations)
                .apply(argsBuilder);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testCalibrationSensitivityThresholdSpecifiedWithoutUnlabeledAnnotations() {
        final File outputDir = createTempDir("train");
        final String outputPrefix = String.format("%s/test", outputDir);
        final ArgumentsBuilder argsBuilder = BASE_ARGUMENTS_BUILDER_SUPPLIER.get();
        argsBuilder.addOutput(outputPrefix);
        final String extractTag = "extract.nonAS.snpIndel.posUn";
        final File positiveAnnotationsHDF5 = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                String.format("%s.annot.hdf5", extractTag));
        final Function<ArgumentsBuilder, ArgumentsBuilder> addPositiveAnnotations = ab ->
                ADD_ANNOTATIONS_HDF5.apply(ab, positiveAnnotationsHDF5);
        final double calibrationSensitivityThreshold = 0.9;
        final Function<ArgumentsBuilder, ArgumentsBuilder> addCalibrationSensitivityThreshold = ab ->
                ADD_CALIBRATION_SENSITIVITY_THRESHOLD.apply(ab, calibrationSensitivityThreshold);
        ADD_ISOLATION_FOREST_PYTHON_SCRIPT
                .andThen(ADD_ISOLATION_FOREST_HYPERPARAMETERS_JSON)
                .andThen(addPositiveAnnotations)
                .andThen(addCalibrationSensitivityThreshold)
                .apply(argsBuilder);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testPositiveAndUnlabeledAnnotationNamesAreNotIdentical() {
        final File outputDir = createTempDir("train");
        final String outputPrefix = String.format("%s/test", outputDir);
        final ArgumentsBuilder argsBuilder = BASE_ARGUMENTS_BUILDER_SUPPLIER.get();
        argsBuilder.addOutput(outputPrefix);
        final File positiveAnnotationsHDF5 = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                "extract.nonAS.snpIndel.posUn.annot.hdf5");         // non-allele-specific
        final Function<ArgumentsBuilder, ArgumentsBuilder> addPositiveAnnotations = ab ->
                ADD_ANNOTATIONS_HDF5.apply(ab, positiveAnnotationsHDF5);
        final File unlabeledAnnotationsHDF5 = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                "extract.AS.snpIndel.posUn.unlabeled.annot.hdf5");  // allele-specific
        final Function<ArgumentsBuilder, ArgumentsBuilder> addUnlabeledAnnotations = ab ->
                ADD_UNLABELED_ANNOTATIONS_HDF5.apply(ab, unlabeledAnnotationsHDF5);
        final double calibrationSensitivityThreshold = 0.9;
        final Function<ArgumentsBuilder, ArgumentsBuilder> addCalibrationSensitivityThreshold = ab ->
                ADD_CALIBRATION_SENSITIVITY_THRESHOLD.apply(ab, calibrationSensitivityThreshold);
        ADD_ISOLATION_FOREST_PYTHON_SCRIPT
                .andThen(ADD_ISOLATION_FOREST_HYPERPARAMETERS_JSON)
                .andThen(addPositiveAnnotations)
                .andThen(addUnlabeledAnnotations)
                .andThen(addCalibrationSensitivityThreshold)
                .apply(argsBuilder);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testPositiveAnnotationsOfSpecifiedVariantTypesNotPresent() {
        final File outputDir = createTempDir("train");
        final String outputPrefix = String.format("%s/test", outputDir);
        final ArgumentsBuilder argsBuilder = BASE_ARGUMENTS_BUILDER_SUPPLIER.get();
        argsBuilder.addOutput(outputPrefix);
        final File positiveAnnotationsHDF5 = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                "extract.nonAS.snp.posUn.annot.hdf5");     // contains only SNPs, but SNP+INDEL is specified
        final Function<ArgumentsBuilder, ArgumentsBuilder> addPositiveAnnotations = ab ->
                ADD_ANNOTATIONS_HDF5.apply(ab, positiveAnnotationsHDF5);
        ADD_SNP_MODE.andThen(ADD_INDEL_MODE)
                .andThen(ADD_ISOLATION_FOREST_PYTHON_SCRIPT)
                .andThen(ADD_ISOLATION_FOREST_HYPERPARAMETERS_JSON)
                .andThen(addPositiveAnnotations)
                .apply(argsBuilder);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testUnlabeledAnnotationsOfSpecifiedVariantTypesNotPresent() {
        final File outputDir = createTempDir("train");
        final String outputPrefix = String.format("%s/test", outputDir);
        final ArgumentsBuilder argsBuilder = BASE_ARGUMENTS_BUILDER_SUPPLIER.get();
        argsBuilder.addOutput(outputPrefix);
        final File positiveAnnotationsHDF5 = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                "extract.nonAS.snpIndel.posUn.annot.hdf5");
        final Function<ArgumentsBuilder, ArgumentsBuilder> addPositiveAnnotations = ab ->
                ADD_ANNOTATIONS_HDF5.apply(ab, positiveAnnotationsHDF5);
        final File unlabeledAnnotationsHDF5 = new File(INPUT_FROM_EXTRACT_EXPECTED_TEST_FILES_DIR,
                "extract.nonAS.snp.posUn.unlabeled.annot.hdf5");    // contains only SNPs, but SNP+INDEL is specified
        final Function<ArgumentsBuilder, ArgumentsBuilder> addUnlabeledAnnotations = ab ->
                ADD_UNLABELED_ANNOTATIONS_HDF5.apply(ab, unlabeledAnnotationsHDF5);
        final double calibrationSensitivityThreshold = 0.9;
        final Function<ArgumentsBuilder, ArgumentsBuilder> addCalibrationSensitivityThreshold = ab ->
                ADD_CALIBRATION_SENSITIVITY_THRESHOLD.apply(ab, calibrationSensitivityThreshold);
        ADD_SNP_MODE.andThen(ADD_INDEL_MODE)
                .andThen(ADD_ISOLATION_FOREST_PYTHON_SCRIPT)
                .andThen(ADD_ISOLATION_FOREST_HYPERPARAMETERS_JSON)
                .andThen(addPositiveAnnotations)
                .andThen(addUnlabeledAnnotations)
                .andThen(addCalibrationSensitivityThreshold)
                .apply(argsBuilder);
        runCommandLine(argsBuilder);
    }
}