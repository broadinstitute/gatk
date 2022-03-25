package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.collect.Lists;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.CommandLineException;
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
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collectors;

/**
 * Note that copies of a subset of the expected outputs for the exact-match tests below are used as inputs for
 * {@link TrainVariantAnnotationsModelIntegrationTest}. Similarly, copies of a subset of the expected outputs for
 * {@link TrainVariantAnnotationsModelIntegrationTest} are used as inputs for {@link ScoreVariantAnnotationsIntegrationTest}.
 * These copies are located in the train/input and score/input subdirectories in the test-resources directory for this package.
 *
 * We choose to use copies of these files (rather than the original expected files) to encapsulate the tests for each tool.
 * However, developers updating any of these tests may optionally choose to keep test files for all tools in sync, as appropriate.
 */
public final class ExtractVariantAnnotationsIntegrationTest extends CommandLineProgramTest {

    // If true, update the expected outputs in tests that assert an exact match vs. prior output,
    // instead of actually running the tests. Can be used with "./gradlew test -Dtest.single=ExtractVariantAnnotationsIntegrationTest"
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

    private static final List<String> NON_ALLELE_SPECIFIC_ANNOTATIONS = Arrays.asList(
            "DP", "FS", "MQ", "MQRankSum", "QD", "ReadPosRankSum", "SOR");

    private static final List<String> ALLELE_SPECIFIC_ANNOTATIONS = Arrays.asList(
            "DP", "AS_FS", "AS_MQ", "AS_MQRankSum", "AS_QD", "AS_ReadPosRankSum", "AS_SOR");

    private static final File TEST_FILES_DIR = new File(largeFileTestDir,
            "org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/");
    private static final File EXPECTED_TEST_FILES_DIR = new File(largeFileTestDir,
            "org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/extract/expected");

    // The input VCF should cover a genomic region given by the union of regions in the below training and truth resources
    // and should also contain a few multiallelics that overlap those resources.
    private static final File INPUT_VCF = new File(TEST_FILES_DIR, "input/stroke_vqsr_magic_as.chr1.1-10M.vcf.gz");

    // We use snippets of the Omni sites for SNP training (chr1:1-5000000) and truth (chr1:5000000-10000000); we don't sweat the 1bp overlap.
    private static final File SNP_TRAINING_VCF = new File(TEST_FILES_DIR, "resources/1000G_omni2.5.hg38.chr1.1-5M.vcf.gz");
    private static final File SNP_TRUTH_VCF = new File(TEST_FILES_DIR, "resources/1000G_omni2.5.hg38.chr1.5M-10M.vcf.gz");

    // We use snippets of the Mills sites for indel training (chr1:1-5000000) and truth (chr1:5000000-10000000); we don't sweat the 1bp overlap.
    private static final File INDEL_TRAINING_VCF = new File(TEST_FILES_DIR, "resources/Mills_and_1000G_gold_standard.indels.hg38.chr1.1-5M.vcf.gz");
    private static final File INDEL_TRUTH_VCF = new File(TEST_FILES_DIR, "resources/Mills_and_1000G_gold_standard.indels.hg38.chr1.5M-10M.vcf.gz");

    private static final int MAXIMUM_NUMBER_OF_UNLABELED_VARIANTS = 100;

    // Supplier and functions for creating and adding various arguments to an ArgumentsBuilder.
    private static final Supplier<ArgumentsBuilder> BASE_ARGS_BUILDER_SUPPLIER = () -> {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.addVCF(INPUT_VCF);
        argsBuilder.add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);
        return argsBuilder;
    };
    private static final Function<ArgumentsBuilder, ArgumentsBuilder> ADD_NON_ALLELE_SPECIFIC_ANNOTATIONS = (argsBuilder) -> {
        NON_ALLELE_SPECIFIC_ANNOTATIONS.forEach(a -> argsBuilder.add(StandardArgumentDefinitions.ANNOTATION_LONG_NAME, a));
        return argsBuilder;
    };
    private static final Function<ArgumentsBuilder, ArgumentsBuilder> ADD_ALLELE_SPECIFIC_ANNOTATIONS = (argsBuilder) -> {
        argsBuilder.addFlag(LabeledVariantAnnotationsWalker.USE_ALLELE_SPECIFIC_ANNOTATIONS_LONG_NAME);
        ALLELE_SPECIFIC_ANNOTATIONS.forEach(a -> argsBuilder.add(StandardArgumentDefinitions.ANNOTATION_LONG_NAME, a));
        return argsBuilder;
    };
    private static final Function<ArgumentsBuilder, ArgumentsBuilder> ADD_SNP_MODE_AND_RESOURCES = (argsBuilder) -> {
        argsBuilder.add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, VariantType.SNP)
                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":omni-training,training=true", SNP_TRAINING_VCF)
                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":omni-truth,truth=true", SNP_TRUTH_VCF);
        return argsBuilder;
    };
    private static final Function<ArgumentsBuilder, ArgumentsBuilder> ADD_INDEL_MODE_AND_RESOURCES = (argsBuilder) -> {
        argsBuilder.add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, VariantType.INDEL)
                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":mills-training,training=true", INDEL_TRAINING_VCF)
                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":mills-truth,truth=true", INDEL_TRUTH_VCF);
        return argsBuilder;
    };
    private static final Function<ArgumentsBuilder, ArgumentsBuilder> ADD_MAXIMUM_NUMBER_OF_UNLABELED_VARIANTS = (argsBuilder) -> {
        argsBuilder.add(ExtractVariantAnnotations.MAXIMUM_NUMBER_OF_UNLABELED_VARIANTS_LONG_NAME, MAXIMUM_NUMBER_OF_UNLABELED_VARIANTS);
        return argsBuilder;
    };

    /**
     * Exact-match tests for all configurations given by the Cartesian product of the following options:
     *  1) non-allele-specific vs. allele-specific
     *  2) SNP vs. indel vs. both
     *  3) positive vs. positive-unlabeled
     */
    @DataProvider(name = "dataValidInputs")
    public Object[][] dataValidInputs() {
        final List<List<Pair<String, Function<ArgumentsBuilder, ArgumentsBuilder>>>> testConfigurations = Lists.cartesianProduct(
                Arrays.asList(
                        Pair.of("nonAS", ADD_NON_ALLELE_SPECIFIC_ANNOTATIONS),
                        Pair.of("AS", ADD_ALLELE_SPECIFIC_ANNOTATIONS)),
                Arrays.asList(
                        Pair.of("snp", ADD_SNP_MODE_AND_RESOURCES),
                        Pair.of("indel", ADD_INDEL_MODE_AND_RESOURCES),
                        Pair.of("both", ADD_SNP_MODE_AND_RESOURCES.andThen(ADD_INDEL_MODE_AND_RESOURCES))),
                Arrays.asList(
                        Pair.of("positive", Function.identity()),
                        Pair.of("positiveUnlabeled", ADD_MAXIMUM_NUMBER_OF_UNLABELED_VARIANTS)));

        return testConfigurations.stream()
                .map(tagAndAddFunctionPairs -> new Object[]{
                        tagAndAddFunctionPairs.stream().map(Pair::getLeft).collect(Collectors.joining(".")), // e.g., nonAS.snp.positive
                        tagAndAddFunctionPairs.stream().map(Pair::getRight)                                              // creates the corresponding ArgumentsBuilder
                                .reduce(Function.identity(), Function::andThen)                                          //  by stringing together functions that add the
                                .apply(BASE_ARGS_BUILDER_SUPPLIER.get())})                                               //  appropriate arguments
                .toArray(Object[][]::new);
    }

    @Test(dataProvider = "dataValidInputs")
    public void testValidInputs(final String tag,
                                final ArgumentsBuilder argsBuilder) {
        final File outputDir = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? EXPECTED_TEST_FILES_DIR : createTempDir("extract");
        final String outputPrefix = String.format("%s/%s", outputDir, tag);
        argsBuilder.addOutput(outputPrefix);
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, "INFO");
        runCommandLine(argsBuilder);

        if (!UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            assertOutputs(tag, outputPrefix);
        }
    }

    private static void assertOutputs(final String tag,
                                      final String outputPrefix) {
        runSystemCommand(String.format("h5diff %s/%s.annot.hdf5 %s.annot.hdf5", EXPECTED_TEST_FILES_DIR, tag, outputPrefix));
        runSystemCommand(String.format("diff %s/%s.vcf.gz %s.vcf.gz", EXPECTED_TEST_FILES_DIR, tag, outputPrefix));
        runSystemCommand(String.format("diff %s/%s.vcf.gz.tbi %s.vcf.gz.tbi", EXPECTED_TEST_FILES_DIR, tag, outputPrefix));
        if (tag.contains("positiveUnlabeled")) {
            runSystemCommand(String.format("h5diff %s/%s.unlabeled.annot.hdf5 %s.unlabeled.annot.hdf5", EXPECTED_TEST_FILES_DIR, tag, outputPrefix));
        } else {
            Assert.assertFalse(new File(outputPrefix, ".unlabeled.annot.hdf5").exists());
        }
    }

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

    @Test(expectedExceptions = CommandLineException.class)
    public void testMissingTrainingResource() {
        final File outputDir = createTempDir("extract");
        final String outputPrefix = String.format("%s/test", outputDir);
        final ArgumentsBuilder argsBuilder = ADD_ALLELE_SPECIFIC_ANNOTATIONS.apply(BASE_ARGS_BUILDER_SUPPLIER.get());
        argsBuilder.add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, "SNP")
                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":omni-training,training=true", SNP_TRAINING_VCF) // only specify training label
                .addOutput(outputPrefix);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testMissingTruthResource() {
        final File outputDir = createTempDir("extract");
        final String outputPrefix = String.format("%s/test", outputDir);
        final ArgumentsBuilder argsBuilder = ADD_ALLELE_SPECIFIC_ANNOTATIONS.apply(BASE_ARGS_BUILDER_SUPPLIER.get());
        argsBuilder.add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, "SNP")
                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":omni-truth,truth=true", SNP_TRUTH_VCF) // only specify truth label
                .addOutput(outputPrefix);
        runCommandLine(argsBuilder);
    }

    // TODO is this expected behavior? or should the tool fail if we don't specify the flag but an AS_* annotation is specified?
    @Test(expectedExceptions = AssertionError.class)
    public void testForgotToSpecifyUseAlleleSpecificAnnotationsFlag() {
        final File outputDir = createTempDir("extract");
        final String outputPrefix = String.format("%s/test", outputDir);
        final ArgumentsBuilder argsBuilder = ADD_SNP_MODE_AND_RESOURCES.apply(BASE_ARGS_BUILDER_SUPPLIER.get());
        ALLELE_SPECIFIC_ANNOTATIONS.forEach(a -> argsBuilder.add(StandardArgumentDefinitions.ANNOTATION_LONG_NAME, a));
        argsBuilder.addOutput(outputPrefix);
        runCommandLine(argsBuilder);

        // check that outputs do not match the expected allele-specific outputs, since we forgot to specify the flag
        assertOutputs("AS.snp.positive", outputPrefix);
    }
}