package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngineArgumentsCollection;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.vcf.VcfUtils;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class SVStratifyIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testBwaMeltCohort() {
        final File outputDir = createTempDir("stratify");
        final String inputVcfPath = getToolTestDataDir() + "bwa_melt.chr22.vcf.gz";
        final String configFile = getToolTestDataDir() + "test_config.tsv";

        final String segdupFile = getToolTestDataDir() + "hg38.SegDup.chr22.bed";
        final String segdupName = "SD";
        final String repeatmaskerFile = getToolTestDataDir() + "hg38.RM.chr22_subsampled.bed";
        final String repeatmaskerName = "RM";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(outputDir)
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test")
                .add(SVStratify.SPLIT_OUTPUT_LONG_NAME, true)
                .add(SVStratificationEngineArgumentsCollection.STRATIFY_CONFIG_FILE_LONG_NAME, configFile)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_NAME_FILE_LONG_NAME, segdupName)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_INTERVAL_FILE_LONG_NAME, segdupFile)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_NAME_FILE_LONG_NAME, repeatmaskerName)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_INTERVAL_FILE_LONG_NAME, repeatmaskerFile)
                .add(SVStratificationEngineArgumentsCollection.OVERLAP_FRACTION_LONG_NAME, 0.5)
                .add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, GATKBaseTest.FULL_HG38_DICT)
                .add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVcfPath);

        runCommandLine(args, SVStratify.class.getSimpleName());

        final File[] outputFiles = outputDir.listFiles();
        Assert.assertEquals(outputFiles.length, 14);
        final Map<String, Integer> expectedOutputSuffixes = new HashMap<>();
        expectedOutputSuffixes.put("INS_small_SD", 46);
        expectedOutputSuffixes.put("DEL_50_5k_both", 110);
        expectedOutputSuffixes.put("DEL_5k_50k_SD", 2);
        expectedOutputSuffixes.put("DUP_lt5kb_RM", 0);
        expectedOutputSuffixes.put("INV_gt1kb", 26);
        expectedOutputSuffixes.put("BND_SD", 77);
        expectedOutputSuffixes.put(SVStratify.DEFAULT_STRATUM, 1196);
        int numVcfs = 0;
        int totalRecords = 0;
        for (final File file : outputFiles) {
            if (VcfUtils.isVariantFile(file)) {
                ++numVcfs;
                final Pair<VCFHeader, List<VariantContext>> outputVcf = VariantContextTestUtils.readEntireVCFIntoMemory(file.getAbsolutePath());
                boolean foundSuffix = false;
                for (final String suffix : expectedOutputSuffixes.keySet()) {
                    if (file.toString().contains("." + suffix + ".")) {
                        foundSuffix = true;
                        for (final VariantContext variant : outputVcf.getRight()) {
                            Assert.assertTrue(variant.hasAttribute(GATKSVVCFConstants.STRATUM_INFO_KEY));
                            Assert.assertEquals(variant.getAttribute(GATKSVVCFConstants.STRATUM_INFO_KEY), suffix);
                        }
                        final int expectedSize = expectedOutputSuffixes.get(suffix).intValue();
                        final int actualSize = outputVcf.getRight().size();
                        Assert.assertEquals(actualSize, expectedSize,
                                "Expected " + expectedSize + " records but found " + actualSize + " in " + suffix);
                        totalRecords += actualSize;
                        break;
                    }
                }
                Assert.assertTrue(foundSuffix, "Unexpected file suffix: " + file.getAbsolutePath());
            }
        }
        Assert.assertEquals(numVcfs, 7);
        final int numInputRecords = VariantContextTestUtils.readEntireVCFIntoMemory(inputVcfPath).getRight().size();
        Assert.assertEquals(totalRecords, numInputRecords);
    }

    @Test
    public void testBwaMeltCohortSingleOutput() {
        final File outputDir = createTempDir("stratify");
        final File outputFile = outputDir.toPath().resolve("out.vcf.gz").toFile();
        final String inputVcfPath = getToolTestDataDir() + "bwa_melt.chr22.vcf.gz";
        final String configFile = getToolTestDataDir() + "test_config.tsv";

        final String segdupFile = getToolTestDataDir() + "hg38.SegDup.chr22.bed";
        final String segdupName = "SD";
        final String repeatmaskerFile = getToolTestDataDir() + "hg38.RM.chr22_subsampled.bed";
        final String repeatmaskerName = "RM";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(outputFile)
                .add(SVStratificationEngineArgumentsCollection.STRATIFY_CONFIG_FILE_LONG_NAME, configFile)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_NAME_FILE_LONG_NAME, segdupName)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_INTERVAL_FILE_LONG_NAME, segdupFile)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_NAME_FILE_LONG_NAME, repeatmaskerName)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_INTERVAL_FILE_LONG_NAME, repeatmaskerFile)
                .add(SVStratificationEngineArgumentsCollection.OVERLAP_FRACTION_LONG_NAME, 0.5)
                .add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, GATKBaseTest.FULL_HG38_DICT)
                .add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVcfPath);

        runCommandLine(args, SVStratify.class.getSimpleName());

        final List<File> outputFiles = Lists.newArrayList(outputDir.listFiles()).stream().filter(VcfUtils::isVariantFile).collect(Collectors.toUnmodifiableList());
        Assert.assertEquals(outputFiles.size(), 1);
        Assert.assertEquals(outputFiles.get(0).getAbsolutePath(), outputFile.getAbsolutePath());
        final Pair<VCFHeader, List<VariantContext>> inputVcf = VariantContextTestUtils.readEntireVCFIntoMemory(inputVcfPath);
        final Pair<VCFHeader, List<VariantContext>> outputVcf = VariantContextTestUtils.readEntireVCFIntoMemory(outputFile.getAbsolutePath());
        Assert.assertEquals(outputVcf.getRight().size(), inputVcf.getRight().size());
    }

    @Test(expectedExceptions = GATKException.class)
    public void testBwaMeltCohortRedundant() {
        final File outputDir = createTempDir("stratify");
        final String inputVcfPath = getToolTestDataDir() + "bwa_melt.chr22.vcf.gz";
        final String configFile = getToolTestDataDir() + "test_config_redundant.tsv";

        final String segdupFile = getToolTestDataDir() + "hg38.SegDup.chr22.bed";
        final String segdupName = "SD";
        final String repeatmaskerFile = getToolTestDataDir() + "hg38.RM.chr22_subsampled.bed";
        final String repeatmaskerName = "RM";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(outputDir)
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test")
                .add(SVStratify.SPLIT_OUTPUT_LONG_NAME, true)
                .add(SVStratificationEngineArgumentsCollection.STRATIFY_CONFIG_FILE_LONG_NAME, configFile)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_NAME_FILE_LONG_NAME, segdupName)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_INTERVAL_FILE_LONG_NAME, segdupFile)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_NAME_FILE_LONG_NAME, repeatmaskerName)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_INTERVAL_FILE_LONG_NAME, repeatmaskerFile)
                .add(SVStratificationEngineArgumentsCollection.OVERLAP_FRACTION_LONG_NAME, 0.5)
                .add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, GATKBaseTest.FULL_HG38_DICT)
                .add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVcfPath);

        runCommandLine(args, SVStratify.class.getSimpleName());
    }

    @Test
    public void testBwaMeltCohortBypassRedundant() {
        final File outputDir = createTempDir("stratify");
        final String inputVcfPath = getToolTestDataDir() + "bwa_melt.chr22.vcf.gz";
        final String configFile = getToolTestDataDir() + "test_config_redundant.tsv";

        final String segdupFile = getToolTestDataDir() + "hg38.SegDup.chr22.bed";
        final String segdupName = "SD";
        final String repeatmaskerFile = getToolTestDataDir() + "hg38.RM.chr22_subsampled.bed";
        final String repeatmaskerName = "RM";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(outputDir)
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test")
                .add(SVStratify.SPLIT_OUTPUT_LONG_NAME, true)
                .add(SVStratificationEngineArgumentsCollection.STRATIFY_CONFIG_FILE_LONG_NAME, configFile)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_NAME_FILE_LONG_NAME, segdupName)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_INTERVAL_FILE_LONG_NAME, segdupFile)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_NAME_FILE_LONG_NAME, repeatmaskerName)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_INTERVAL_FILE_LONG_NAME, repeatmaskerFile)
                .add(SVStratificationEngineArgumentsCollection.OVERLAP_FRACTION_LONG_NAME, 0.5)
                .add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, GATKBaseTest.FULL_HG38_DICT)
                .add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVcfPath)
                .addFlag(SVStratify.ALLOW_MULTIPLE_MATCHES_LONG_NAME);

        runCommandLine(args, SVStratify.class.getSimpleName());
    }

    @Test(expectedExceptions = {GATKException.class, IllegalArgumentException.class})
    public void testBwaMeltCohortDuplicateContextName() {
        final File outputDir = createTempDir("stratify");
        final String inputVcfPath = getToolTestDataDir() + "bwa_melt.chr22.vcf.gz";
        final String configFile = getToolTestDataDir() + "test_config_duplicate.tsv";

        final String segdupFile = getToolTestDataDir() + "hg38.SegDup.chr22.bed";
        final String segdupName = "SD";
        final String repeatmaskerFile = getToolTestDataDir() + "hg38.RM.chr22_subsampled.bed";
        final String repeatmaskerName = "RM";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(outputDir)
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, "test")
                .add(SVStratify.SPLIT_OUTPUT_LONG_NAME, true)
                .add(SVStratificationEngineArgumentsCollection.STRATIFY_CONFIG_FILE_LONG_NAME, configFile)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_NAME_FILE_LONG_NAME, segdupName)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_INTERVAL_FILE_LONG_NAME, segdupFile)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_NAME_FILE_LONG_NAME, repeatmaskerName)
                .add(SVStratificationEngineArgumentsCollection.CONTEXT_INTERVAL_FILE_LONG_NAME, repeatmaskerFile)
                .add(SVStratificationEngineArgumentsCollection.OVERLAP_FRACTION_LONG_NAME, 0.5)
                .add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, GATKBaseTest.FULL_HG38_DICT)
                .add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVcfPath);

        runCommandLine(args, SVStratify.class.getSimpleName());
    }
}