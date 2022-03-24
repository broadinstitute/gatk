package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.collect.Lists;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collectors;

/**
 * TODO
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
            "org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/expected/extract");
    private static final File INPUT_VCF = new File(TEST_FILES_DIR, "input/stroke_vqsr_magic_as.chr1.1-10M.vcf.gz");
    private static final File SNP_TRAINING_VCF = new File(TEST_FILES_DIR, "resources/1000G_omni2.5.hg38.chr1.1-5M.vcf.gz");
    private static final File SNP_TRUTH_VCF = new File(TEST_FILES_DIR, "resources/1000G_omni2.5.hg38.chr1.5M-10M.vcf.gz");
    private static final File INDEL_TRAINING_VCF = new File(TEST_FILES_DIR, "resources/Mills_and_1000G_gold_standard.indels.hg38.chr1.1-5M.vcf.gz");
    private static final File INDEL_TRUTH_VCF = new File(TEST_FILES_DIR, "resources/Mills_and_1000G_gold_standard.indels.hg38.chr1.5M-10M.vcf.gz");
    private static final int MAXIMUM_NUMBER_OF_UNLABELED_VARIANTS = 100;

    @DataProvider(name = "dataValidInputs")
    public Object[][] dataValidInputs() {
        final Supplier<ArgumentsBuilder> baseArgsBuilderSupplier = () -> {
            final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
            argsBuilder.addVCF(INPUT_VCF);
            return argsBuilder;
        };

        final Function<ArgumentsBuilder, ArgumentsBuilder> addNonAlleleSpecificAnnotations = (argsBuilder) -> {
            NON_ALLELE_SPECIFIC_ANNOTATIONS.forEach(a -> argsBuilder.add(StandardArgumentDefinitions.ANNOTATION_LONG_NAME, a));
            return argsBuilder;
        };
        final Function<ArgumentsBuilder, ArgumentsBuilder> addAlleleSpecificAnnotations = (argsBuilder) -> {
            argsBuilder.addFlag(LabeledVariantAnnotationsWalker.USE_ALLELE_SPECIFIC_ANNOTATIONS_LONG_NAME); // TODO check why removing this line still works
            ALLELE_SPECIFIC_ANNOTATIONS.forEach(a -> argsBuilder.add(StandardArgumentDefinitions.ANNOTATION_LONG_NAME, a));
            return argsBuilder;
        };

        final Function<ArgumentsBuilder, ArgumentsBuilder> addSNPModeAndResources = (argsBuilder) -> {
            argsBuilder.add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, "SNP")
                    .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":omni-training,training=true", SNP_TRAINING_VCF)
                    .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":omni-truth,truth=true", SNP_TRUTH_VCF);
            return argsBuilder;
        };
        final Function<ArgumentsBuilder, ArgumentsBuilder> addIndelModeAndResources = (argsBuilder) -> {
            argsBuilder.add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, "INDEL")
                    .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":mills-training,training=true", INDEL_TRAINING_VCF)
                    .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":mills-truth,truth=true", INDEL_TRUTH_VCF);
            return argsBuilder;
        };

        final Function<ArgumentsBuilder, ArgumentsBuilder> addMaximumNumberOfUnlabeledVariants = (argsBuilder) -> {
            argsBuilder.add(ExtractVariantAnnotations.MAXIMUM_NUMBER_OF_UNLABELED_VARIANTS_LONG_NAME, MAXIMUM_NUMBER_OF_UNLABELED_VARIANTS);
            return argsBuilder;
        };

        final List<List<Pair<String, Function<ArgumentsBuilder, ArgumentsBuilder>>>> testConfigurations = Lists.cartesianProduct(
                Arrays.asList(
                        Pair.of("nonAS", addNonAlleleSpecificAnnotations),
                        Pair.of("AS", addAlleleSpecificAnnotations)),
                Arrays.asList(
                        Pair.of("snp", addSNPModeAndResources),
                        Pair.of("indel", addIndelModeAndResources),
                        Pair.of("both", addSNPModeAndResources.andThen(addIndelModeAndResources))),
                Arrays.asList(
                        Pair.of("positive", Function.identity()),
                        Pair.of("positiveUnlabeled", addMaximumNumberOfUnlabeledVariants)));

        return testConfigurations.stream()
                .map(tagAndAddFunctionPairs -> new Object[]{
                        tagAndAddFunctionPairs.stream().map(Pair::getLeft).collect(Collectors.joining(".")),
                        tagAndAddFunctionPairs.stream().map(Pair::getRight)
                                .reduce(Function.identity(), Function::andThen)
                                .apply(baseArgsBuilderSupplier.get())})
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

    }

//    private static void assertAnnotationFilesEqual(final File annotationsFile1,
//                                                   final File annotationsFile2) {
//        try (final HDF5File annotationsHDF5File1 = new HDF5File(annotationsFile1, HDF5File.OpenMode.READ_ONLY);
//             final HDF5File annotationsHDF5File2 = new HDF5File(annotationsFile2, HDF5File.OpenMode.READ_ONLY)) {
//            Assert.assertEquals(
//                    annotationsHDF5File1.readStringArray("/annotations/names"),
//                    annotationsHDF5File2.readStringArray("/annotations/names"));
//            for (final String label : Arrays.asList("training", "truth")) {
//                Assert.assertEquals(
//                        LabeledVariantAnnotationsData.readChunkedDoubleArray(annotationsHDF5File1, String.format("/labels/%s", label)),
//                        LabeledVariantAnnotationsData.readChunkedDoubleArray(annotationsHDF5File2, String.format("/labels/%s", label)));
//            }
//            Assert.assertEquals(
//                    HDF5Utils.readChunkedDoubleMatrix(annotationsHDF5File1, "/annotations"),
//                    HDF5Utils.readChunkedDoubleMatrix(annotationsHDF5File2, "/annotations"));
//        } catch (final HDF5LibException exception) {
//            Assert.fail("Exception encountered during reading of annotations:", exception);
//        }
//    }

    @Test
    public void test1kgp50ExomesAllUnlabeled() throws IOException {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/extract-test/test.all-unlabeled",
                "--maximum-number-of-unlabeled-variants", "10000",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "--mode", "SNP",
                "--mode", "INDEL",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--resource:mills,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);

        runSystemCommand("h5diff /home/slee/working/vqsr/scalable/extract-test/test.all-unlabeled.annot.hdf5 /home/slee/working/vqsr/scalable/extract-test/expected/test.all-unlabeled.annot.hdf5");
        runSystemCommand("h5diff /home/slee/working/vqsr/scalable/extract-test/test.all-unlabeled.unlabeled.annot.hdf5 /home/slee/working/vqsr/scalable/extract-test/expected/test.all-unlabeled.unlabeled.annot.hdf5");
        runSystemCommand("diff /home/slee/working/vqsr/scalable/extract-test/test.all-unlabeled.vcf.gz /home/slee/working/vqsr/scalable/extract-test/expected/test.all-unlabeled.vcf.gz");
        runSystemCommand("diff /home/slee/working/vqsr/scalable/extract-test/test.all-unlabeled.vcf.gz.tbi /home/slee/working/vqsr/scalable/extract-test/expected/test.all-unlabeled.vcf.gz.tbi");
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
            e.printStackTrace();
        }
    }

    @Test
    public void test1kgp50ExomesAll() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/extract-test/test.all",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "--trust-all-polymorphic",
                "--mode", "SNP",
                "--mode", "INDEL",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--resource:mills,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesSNP() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/extract-test/test.snp",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "--trust-all-polymorphic",
                "--mode", "SNP",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testSNPAS() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", largeFileTestDir + "VQSR/chr1snippet.doctoredMQ.doctoredAS.sites_only.vcf",
                "-O", "/home/slee/working/vqsr/scalable/extract-test/test.snp.as",
                "--use-allele-specific-annotations",
                "-A", "AS_FS",
                "-A", "AS_ReadPosRankSum",
                "-A", "AS_MQRankSum",
                "-A", "AS_QD",
                "-A", "AS_SOR",
                "-A", "AS_MQ",
                "--trust-all-polymorphic",
                "--mode", "SNP",
                "--resource:doctored,training=true,truth=true", largeFileTestDir + "VQSR/chr1snippet.doctoredMQ.doctoredAS.sites_only.vcf",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testSNPNonAS() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", largeFileTestDir + "VQSR/chr1snippet.doctoredMQ.doctoredAS.sites_only.vcf",
                "-O", "/home/slee/working/vqsr/scalable/extract-test/test.snp.non-as",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "--trust-all-polymorphic",
                "--mode", "SNP",
                "--resource:doctored,training=true,truth=true", largeFileTestDir + "VQSR/chr1snippet.doctoredMQ.doctoredAS.sites_only.vcf",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesIndel() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/extract-test/test.indel",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "--trust-all-polymorphic",
                "--mode", "INDEL",
                "--resource:mills,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testJbxAll() {
        final String[] arguments = {
                "-V", "/home/slee/working/vqsr/scalable/jbx/resources/Test50Callset.annoated_pids.sites-only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all.extract",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "-A", "COMBINED_TREE_SCORE",
                "--trust-all-polymorphic",
                "--mode", "SNP",
                "--mode", "INDEL",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--resource:mills,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testJbxAllUnlabeled() {
        final String[] arguments = {
                "-V", "/home/slee/working/vqsr/scalable/jbx/resources/Test50Callset.annoated_pids.sites-only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all-unlabeled.extract",
                "--maximum-number-of-unlabeled-variants", "10000000",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "-A", "COMBINED_TREE_SCORE",
                "--trust-all-polymorphic",
                "--mode", "SNP",
                "--mode", "INDEL",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--resource:mills,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testJbxSNP() {
        final String[] arguments = {
                "-V", "/home/slee/working/vqsr/scalable/jbx/resources/Test50Callset.annoated_pids.sites-only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.snp.extract",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "-A", "COMBINED_TREE_SCORE",
                "--trust-all-polymorphic",
                "--mode", "SNP",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }
}