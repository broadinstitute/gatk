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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
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
    private static final File EXPECTED_TEST_FILES_DIR = new File(TEST_FILES_DIR, "expected");

    private static final File ISOLATION_FOREST_PYTHON_SCRIPT = new File(packageMainResourcesDir,
            "tools/walkers/vqsr/scalable/isolation-forest.py");
    private static final File ISOLATION_FOREST_HYPERPARAMETERS_JSON = new File(TEST_FILES_DIR,
            "input/isolation-forest-hyperparameters.json");

    // Supplier and functions for creating and adding various arguments to an ArgumentsBuilder.
    private static final Supplier<ArgumentsBuilder> BASE_ARGS_BUILDER_SUPPLIER = ArgumentsBuilder::new;

    private static final BiFunction<ArgumentsBuilder, File, ArgumentsBuilder> ADD_ANNOTATIONS_HDF5 = (argsBuilder, annotationsHDF5) -> {
        argsBuilder.add(TrainVariantAnnotationsModel.ANNOTATIONS_HDF5_LONG_NAME, annotationsHDF5);
        return argsBuilder;
    };
    private static final BiFunction<ArgumentsBuilder, File, ArgumentsBuilder> ADD_UNLABELED_ANNOTATIONS_HDF5 = (argsBuilder, unlabeledAnnotationsHDF5) -> {
        argsBuilder.add(TrainVariantAnnotationsModel.UNLABELED_ANNOTATIONS_HDF5_LONG_NAME, unlabeledAnnotationsHDF5);
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
     * Exact-match tests for all configurations given by the Cartesian product of the following options:
     *  1) SNP-only vs. SNP+INDEL (for both of these options, we use extracted annotations that contain both SNP and INDEL variants as input)
     *  2) positive (training with nonAS.both.positive.annot.hdf5) vs. positive-unlabeled (training with both nonAS.both.positive.annot.hdf5 and nonAS.both.unlabeled.annot.hdf5)
     *  3) Java Bayesian Gaussian Mixture Model (BGMM) backend vs. python sklearn IsolationForest backend TODO the BGMM has been reduced to a stub for this initial PR; subsequent PRs will cover the backend code and reconnect the stub
     */
    @DataProvider(name = "dataValidInputs")
    public Object[][] dataValidInputs() {
        final File positiveAnnotationsHDF5 = new File(TEST_FILES_DIR, "input/extract.nonAS.snpIndel.posUn.annot.hdf5");
        final File unlabeledAnnotationsHDF5 = new File(TEST_FILES_DIR, "input/extract.nonAS.snpIndel.posUn.unlabeled.annot.hdf5");
        final Function<ArgumentsBuilder, ArgumentsBuilder> addPositiveAnnotations = argsBuilder ->
                ADD_ANNOTATIONS_HDF5.apply(argsBuilder, positiveAnnotationsHDF5);
        final Function<ArgumentsBuilder, ArgumentsBuilder> addUnlabeledAnnotations = argsBuilder ->
                ADD_UNLABELED_ANNOTATIONS_HDF5.apply(argsBuilder, unlabeledAnnotationsHDF5);
        final List<List<Pair<String, Function<ArgumentsBuilder, ArgumentsBuilder>>>> testConfigurations = Lists.cartesianProduct(
                Arrays.asList(
                        Pair.of("extract.nonAS.snpIndel.posUn.train.snp", ADD_SNP_MODE),
                        Pair.of("extract.nonAS.snpIndel.posUn.train.snpIndel", ADD_SNP_MODE.andThen(ADD_INDEL_MODE))),
                Arrays.asList(
                        Pair.of("pos", addPositiveAnnotations),
                        Pair.of("posNeg", addPositiveAnnotations.andThen(addUnlabeledAnnotations))),
                Arrays.asList(
//                        Pair.of("BGMM", Function.identity()),
                        Pair.of("IF", ADD_ISOLATION_FOREST_PYTHON_SCRIPT.andThen(ADD_ISOLATION_FOREST_HYPERPARAMETERS_JSON))));

        return testConfigurations.stream()
                .map(tagAndAddFunctionPairs -> new Object[]{
                        tagAndAddFunctionPairs.stream().map(Pair::getLeft).collect(Collectors.joining(".")), // e.g., extract.nonAS.snpIndel.posUn.train.snp.pos.IF
                        tagAndAddFunctionPairs.stream().map(Pair::getRight)                                              // creates the corresponding ArgumentsBuilder
                                .reduce(Function.identity(), Function::andThen)                                          //  by stringing together functions that add the
                                .apply(BASE_ARGS_BUILDER_SUPPLIER.get())})                                               //  appropriate arguments
                .toArray(Object[][]::new);
    }

    @Test(dataProvider = "dataValidInputs")
    public void testValidInputs(final String tag,
                                final ArgumentsBuilder argsBuilder) {
        final File outputDir = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? EXPECTED_TEST_FILES_DIR : createTempDir("train");
        final String outputPrefix = String.format("%s/%s", outputDir, tag);
        argsBuilder.addOutput(outputPrefix);
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, "INFO");
        runCommandLine(argsBuilder);

        if (!UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            assertOutputs(tag, outputPrefix);
        }
    }

    // test using nonAS.snp.positive for SNP-only training = using nonAS.both.positive for SNP-only training

    // exception test annotation-name validation using AS.snp.positive and nonAS.snp.positiveUnlabeled
    // exception test no-data validation using nonAS.snp.positive and nonAS.both.positiveUnlabeled for SNP+INDEL training
    // exception test no-data validation using nonAS.snp.positive for INDEL training

    private static void assertOutputs(final String tag,
                                      final String outputPrefix) {
        if (tag.contains("train.snp")) {
            assertVariantTypeOutputs(tag, outputPrefix, "snp");
            assertNotVariantTypeOutputs(outputPrefix, "indel");
        } else if (tag.contains("train.snpIndel")) {
            assertVariantTypeOutputs(tag, outputPrefix, "snp");
            assertVariantTypeOutputs(tag, outputPrefix, "indel");
        } else {
            Assert.fail("Unknown variant-type tag.");
        }
    }

    private static void assertVariantTypeOutputs(final String tag,
                                                 final String outputPrefix,
                                                 final String variantType) {
        final String tagAndVariantType = String.format("%s.%s", tag, variantType);
        final String outputPrefixAndVariantType = String.format("%s.%s", outputPrefix, variantType);
        runSystemCommand(String.format("h5diff %s/%s.trainingScores.hdf5 %s.trainingScores.hdf5", EXPECTED_TEST_FILES_DIR, tagAndVariantType, outputPrefixAndVariantType));
        runSystemCommand(String.format("h5diff %s/%s.truthScores.hdf5 %s.truthScores.hdf5", EXPECTED_TEST_FILES_DIR, tagAndVariantType, outputPrefixAndVariantType));
        if (tag.contains("BGMM")) {
            runSystemCommand(String.format("diff %s/%s.scorer.ser %s.scorer.ser", EXPECTED_TEST_FILES_DIR, tagAndVariantType, outputPrefixAndVariantType));
            Assert.assertFalse(new File(outputPrefixAndVariantType, ".scorer.pkl").exists());
        } else if (tag.contains("IF")) {
            runSystemCommand(String.format("diff %s/%s.scorer.pkl %s.scorer.pkl", EXPECTED_TEST_FILES_DIR, tagAndVariantType, outputPrefixAndVariantType));
            Assert.assertFalse(new File(outputPrefixAndVariantType, ".scorer.ser").exists());
        } else {
            Assert.fail("Unknown model-backend tag.");
        }
    }

    private static void assertNotVariantTypeOutputs(final String outputPrefix,
                                                    final String variantType) {
        final String outputPrefixAndVariantType = String.format("%s.%s", outputPrefix, variantType);
        Assert.assertFalse(new File(outputPrefixAndVariantType, ".trainingScores.hdf5").exists());
        Assert.assertFalse(new File(outputPrefixAndVariantType, ".truthScores.hdf5").exists());
        Assert.assertFalse(new File(outputPrefixAndVariantType, ".scorer.ser").exists());
        Assert.assertFalse(new File(outputPrefixAndVariantType, ".scorer.pkl").exists());
    }

    // this method is duplicated in the other integration-test classes in this package
    private static void runSystemCommand(final String command) {
        try {
            final Process process = Runtime.getRuntime().exec(command);
            final BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            while (reader.readLine() != null) {
                Assert.fail(command);
            }
            reader.close();
        } catch (final IOException e) {
            Assert.fail(e.getMessage());
        }
    }

//    @Test
//    public void test1kgp50ExomesAll() {
//        final String[] arguments = {
//                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.all.annot.hdf5",
//                "-O", "/home/slee/working/vqsr/scalable/train-test/test.all",
//                "--python-script", PYTHON_SCRIPT,
//                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/hyperparameters.json",
//                "--mode", "SNP",
//                "--mode", "INDEL",
//                "--verbosity", "INFO"
//        };
//        runCommandLine(arguments);
//    }
//
//    @Test
//    public void test1kgp50ExomesAllUnlabeled() {
//        final String[] arguments = {
//                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.all-unlabeled.annot.hdf5",
//                "--unlabeled-annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.all-unlabeled.unlabeled.annot.hdf5",
//                "--truth-sensitivity-threshold", "0.95",
//                "-O", "/home/slee/working/vqsr/scalable/train-test/test.all-unlabeled",
//                "--python-script", PYTHON_SCRIPT,
//                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/hyperparameters.json",
//                "--mode", "SNP",
//                "--mode", "INDEL",
//                "--verbosity", "INFO"
//        };
//        runCommandLine(arguments);
//
//        runSystemCommand("h5diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.snp.trainingScores.hdf5 /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.snp.trainingScores.hdf5");
//        runSystemCommand("h5diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.snp.truthScores.hdf5 /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.snp.truthScores.hdf5");
//        runSystemCommand("h5diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.snp.unlabeledScores.hdf5 /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.snp.unlabeledScores.hdf5");
//        runSystemCommand("diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.snp.scorer.pkl /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.snp.scorer.pkl");
//        runSystemCommand("diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.snp.negative.scorer.pkl /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.snp.negative.scorer.pkl");
//        runSystemCommand("h5diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.indel.trainingScores.hdf5 /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.indel.trainingScores.hdf5");
//        runSystemCommand("h5diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.indel.truthScores.hdf5 /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.indel.truthScores.hdf5");
//        runSystemCommand("h5diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.indel.unlabeledScores.hdf5 /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.indel.unlabeledScores.hdf5");
//        runSystemCommand("diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.indel.scorer.pkl /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.indel.scorer.pkl");
//        runSystemCommand("diff /home/slee/working/vqsr/scalable/train-test/test.all-unlabeled.indel.negative.scorer.pkl /home/slee/working/vqsr/scalable/train-test/expected/test.all-unlabeled.indel.negative.scorer.pkl");
//    }
//
//    @Test
//    public void test1kgp50ExomesSNP() {
//        final String[] arguments = {
//                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.snp.annot.hdf5",
//                "-O", "/home/slee/working/vqsr/scalable/train-test/test",
//                "--python-script", PYTHON_SCRIPT,
//                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/hyperparameters.json",
//                "--mode", "SNP",
//                "--verbosity", "INFO"
//        };
//        runCommandLine(arguments);
//    }
//
//    @Test
//    public void test1kgp50ExomesIndel() {
//        final String[] arguments = {
//                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.indel.annot.hdf5",
//                "-O", "/home/slee/working/vqsr/scalable/train-test/test",
//                "--python-script", PYTHON_SCRIPT,
//                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/hyperparameters.json",
//                "--mode", "INDEL",
//                "--verbosity", "INFO"
//        };
//        runCommandLine(arguments);
//    }
//
//    @Test
//    public void testSNPAS() {
//        final String[] arguments = {
//                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.snp.as.annot.hdf5",
//                "-O", "/home/slee/working/vqsr/scalable/train-test/test.snp.as",
//                "--python-script", PYTHON_SCRIPT,
//                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/hyperparameters.json",
//                "--mode", "SNP",
//                "--verbosity", "INFO"
//        };
//        runCommandLine(arguments);
//    }
//
//    @Test
//    public void test1kgp50ExomesBGMMAll() {
//        final String[] arguments = {
//                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.all.annot.hdf5",
//                "-O", "/home/slee/working/vqsr/scalable/train-test/test.bgmm.all",
//                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/bgmm-hyperparameters.json",
//                "--mode", "SNP",
//                "--mode", "INDEL",
//                "--verbosity", "INFO"
//        };
//        runCommandLine(arguments);
//    }
//
//    @Test
//    public void test1kgp50ExomesBGMMSNP() {
//        final String[] arguments = {
//                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.snp.annot.hdf5",
//                "-O", "/home/slee/working/vqsr/scalable/train-test/test.bgmm",
//                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/bgmm-hyperparameters.json",
//                "--mode", "SNP",
//                "--verbosity", "INFO"
//        };
//        runCommandLine(arguments);
//    }
//
//    @Test
//    public void test1kgp50ExomesBGMMIndel() {
//        final String[] arguments = {
//                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/extract-test/test.indel.annot.hdf5",
//                "-O", "/home/slee/working/vqsr/scalable/train-test/test.bgmm",
//                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/bgmm-hyperparameters.json",
//                "--mode", "INDEL",
//                "--verbosity", "INFO"
//        };
//        runCommandLine(arguments);
//    }
//
//    @Test
//    public void testJbxAll() {
//        final String[] arguments = {
//                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all.extract.annot.hdf5",
//                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all.train",
//                "--python-script", PYTHON_SCRIPT,
//                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/jbx/hyperparameters.json",
//                "--mode", "SNP",
//                "--mode", "INDEL",
//                "--verbosity", "INFO"
//        };
//        runCommandLine(arguments);
//    }
//
//    @Test
//    public void testJbxAllUnlabeled() {
//        final String[] arguments = {
//                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all-unlabeled.extract.annot.hdf5",
//                "--unlabeled-annotations-hdf5", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all-unlabeled.extract.unlabeled.annot.hdf5",
//                "--truth-sensitivity-threshold", "0.95",
//                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all-unlabeled.train",
//                "--python-script", PYTHON_SCRIPT,
//                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/jbx/hyperparameters.json",
//                "--mode", "SNP",
//                "--mode", "INDEL",
//                "--verbosity", "INFO"
//        };
//        runCommandLine(arguments);
//    }
//
//    @Test
//    public void testJbxSNP() {
//        final String[] arguments = {
//                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.snp.extract.annot.hdf5",
//                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.snp.train",
//                "--python-script", PYTHON_SCRIPT,
//                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/jbx/hyperparameters.json",
//                "--mode", "SNP",
//                "--verbosity", "INFO"
//        };
//        runCommandLine(arguments);
//    }
//
//    @Test
//    public void testJbxBGMMAll() {
//        final String[] arguments = {
//                "--annotations-hdf5", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all.extract.annot.hdf5",
//                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.bgmm.all.train",
//                "--hyperparameters-json", "/home/slee/working/vqsr/scalable/train-test/bgmm-hyperparameters.json",
//                "--mode", "SNP",
//                "--mode", "INDEL",
//                "--verbosity", "INFO"
//        };
//        runCommandLine(arguments);
//    }
}