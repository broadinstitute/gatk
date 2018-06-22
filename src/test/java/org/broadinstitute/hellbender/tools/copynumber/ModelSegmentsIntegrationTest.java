package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.*;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.models.AlleleFractionParameter;
import org.broadinstitute.hellbender.tools.copynumber.models.CopyRatioParameter;
import org.broadinstitute.hellbender.tools.copynumber.models.MultidimensionalModellerUnitTest;
import org.broadinstitute.hellbender.tools.copynumber.segmentation.MultidimensionalKernelSegmenterUnitTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

/**
 * Integration tests for {@link ModelSegments}.  We test for input validation across various run modes of the tool
 * and for consistency of metadata in output, but do not test for correctness of output (which is tested elsewhere,
 * e.g. {@link MultidimensionalKernelSegmenterUnitTest} and {@link MultidimensionalModellerUnitTest}).
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class ModelSegmentsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber");
    private static final File TUMOR_DENOISED_COPY_RATIOS_FILE = new File(TEST_SUB_DIR,
            "model-segments-wes-tumor-denoised-copy-ratios-SM-74P4M-v1-chr20-downsampled.deduplicated.denoisedCR.tsv");
    private static final File TUMOR_ALLELIC_COUNTS_FILE = new File(TEST_SUB_DIR,
            "model-segments-wes-tumor-allelic-counts-SM-74P4M-v1-chr20-downsampled.deduplicated.allelicCounts.tsv");
    private static final File NORMAL_ALLELIC_COUNTS_FILE = new File(TEST_SUB_DIR,
            "model-segments-wes-normal-allelic-counts-SM-74NEG-v1-chr20-downsampled.deduplicated.allelicCounts.tsv");
    private static final File TUMOR_DENOISED_COPY_RATIOS_WITH_SAMPLE_NAME_MISMATCH_FILE = new File(TEST_SUB_DIR,
            "model-segments-wes-tumor-denoised-copy-ratios-with-sample-name-mismatch.denoisedCR.tsv");
    private static final File NORMAL_ALLELIC_COUNTS_FILE_WITH_MISSING_SITES = new File(TEST_SUB_DIR,
            "model-segments-wes-normal-allelic-counts-with-missing-sites.allelicCounts.tsv");

    private static final SampleLocatableMetadata EXPECTED_METADATA = new CopyRatioCollection(TUMOR_DENOISED_COPY_RATIOS_FILE).getMetadata();

    @Test
    public void testAllInputsAvailable() {
        final File outputDir = createTempDir("testDir");
        final String outputPrefix = "test";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, TUMOR_DENOISED_COPY_RATIOS_FILE.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME, NORMAL_ALLELIC_COUNTS_FILE.getAbsolutePath())
                .addOutput(outputDir)
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, outputPrefix);
        runCommandLine(argsBuilder);
        assertOutputFiles(outputDir, outputPrefix, true, true);
    }

    @Test
    public void testNoNormalAllelicCounts() {
        final File outputDir = createTempDir("testDir");
        final String outputPrefix = "test";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, TUMOR_DENOISED_COPY_RATIOS_FILE.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath())
                .addOutput(outputDir)
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, outputPrefix);
        runCommandLine(argsBuilder);
        assertOutputFiles(outputDir, outputPrefix, true, false);
    }

    @Test
    public void testNoDenoisedCopyRatios() {
        final File outputDir = createTempDir("testDir");
        final String outputPrefix = "test";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME, NORMAL_ALLELIC_COUNTS_FILE.getAbsolutePath())
                .addOutput(outputDir)
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, outputPrefix);
        runCommandLine(argsBuilder);
        assertOutputFiles(outputDir, outputPrefix, true, true);
    }

    @Test
    public void testAllelicCountsOnly() {
        final File outputDir = createTempDir("testDir");
        final String outputPrefix = "test";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath())
                .addOutput(outputDir)
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, outputPrefix);
        runCommandLine(argsBuilder);
        assertOutputFiles(outputDir, outputPrefix, true, false);
    }

    @Test
    public void testNoAllelicCounts() {
        final File outputDir = createTempDir("testDir");
        final String outputPrefix = "test";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, TUMOR_DENOISED_COPY_RATIOS_FILE.getAbsolutePath())
                .addOutput(outputDir)
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, outputPrefix);
        runCommandLine(argsBuilder);
        assertOutputFiles(outputDir, outputPrefix, false, false);
    }

    private static void assertOutputFiles(final File outputDir,
                                          final String outputPrefix,
                                          final boolean isAllelicCountsPresent,
                                          final boolean isNormalAllelicCountsPresent) {
        Assert.assertTrue(!(!isAllelicCountsPresent && isNormalAllelicCountsPresent));
        for (final String fileTag : Arrays.asList(ModelSegments.BEGIN_FIT_FILE_TAG, ModelSegments.FINAL_FIT_FILE_TAG)) {
            final ModeledSegmentCollection modeledSegments =
                    new ModeledSegmentCollection(new File(outputDir, outputPrefix + fileTag + ModelSegments.SEGMENTS_FILE_SUFFIX));
            Assert.assertEquals(EXPECTED_METADATA, modeledSegments.getMetadata());

            final ParameterDecileCollection<CopyRatioParameter> copyRatioParameters = new ParameterDecileCollection<>(new File(outputDir,
                    outputPrefix + fileTag + ModelSegments.COPY_RATIO_MODEL_PARAMETER_FILE_SUFFIX), CopyRatioParameter.class);
            Assert.assertEquals(EXPECTED_METADATA.getSampleName(), copyRatioParameters.getMetadata().getSampleName());

            final ParameterDecileCollection<AlleleFractionParameter> alleleFractionParameters = new ParameterDecileCollection<>(new File(outputDir,
                    outputPrefix + fileTag + ModelSegments.ALLELE_FRACTION_MODEL_PARAMETER_FILE_SUFFIX), AlleleFractionParameter.class);
            Assert.assertEquals(EXPECTED_METADATA.getSampleName(), alleleFractionParameters.getMetadata().getSampleName());
        }

        final CopyRatioSegmentCollection copyRatioSegments = new CopyRatioSegmentCollection(new File(outputDir,
                outputPrefix + ModelSegments.COPY_RATIO_SEGMENTS_FOR_CALLER_FILE_SUFFIX));
        Assert.assertEquals(EXPECTED_METADATA, copyRatioSegments.getMetadata());

        final File copyRatioLegacySegmentsFile = new File(outputDir,
                outputPrefix + ModelSegments.COPY_RATIO_LEGACY_SEGMENTS_FILE_SUFFIX);
        Assert.assertTrue(copyRatioLegacySegmentsFile.exists());
        final File alleleFractionLegacySegmentsFile = new File(outputDir,
                outputPrefix + ModelSegments.ALLELE_FRACTION_LEGACY_SEGMENTS_FILE_SUFFIX);
        Assert.assertTrue(alleleFractionLegacySegmentsFile.exists());

        AllelicCountCollection hetAllelicCounts = null;
        AllelicCountCollection hetNormalAllelicCounts = null;
        if (isAllelicCountsPresent) {
            hetAllelicCounts = new AllelicCountCollection(new File(outputDir,
                    outputPrefix + ModelSegments.HET_ALLELIC_COUNTS_FILE_SUFFIX));
            Assert.assertEquals(EXPECTED_METADATA, hetAllelicCounts.getMetadata());
        }
        if (isNormalAllelicCountsPresent) { //if this is true, case sample allelic counts will be present
            hetNormalAllelicCounts = new AllelicCountCollection(new File(outputDir,
                    outputPrefix + ModelSegments.NORMAL_HET_ALLELIC_COUNTS_FILE_SUFFIX));
            Assert.assertNotEquals(EXPECTED_METADATA, hetNormalAllelicCounts.getMetadata());    //sample names should differ
            Assert.assertTrue(CopyNumberArgumentValidationUtils.isSameDictionary(               //sequence dictionary should be the same
                    EXPECTED_METADATA.getSequenceDictionary(), hetNormalAllelicCounts.getMetadata().getSequenceDictionary()));
            Assert.assertEquals(hetAllelicCounts.getIntervals(), hetNormalAllelicCounts.getIntervals());
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testSampleNameMismatch() {
        final File outputDir = createTempDir("testDir");
        final String outputPrefix = "test";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, TUMOR_DENOISED_COPY_RATIOS_WITH_SAMPLE_NAME_MISMATCH_FILE.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME, NORMAL_ALLELIC_COUNTS_FILE.getAbsolutePath())
                .addOutput(outputDir)
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, outputPrefix);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = UserException.class)
    public void testMissingSites() {
        final File outputDir = createTempDir("testDir");
        final String outputPrefix = "test";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME, NORMAL_ALLELIC_COUNTS_FILE_WITH_MISSING_SITES.getAbsolutePath())
                .addOutput(outputDir)
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, outputPrefix);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = UserException.class)
    public void testOutputDirExists() {
        final String outputPrefix = "test";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, TUMOR_DENOISED_COPY_RATIOS_FILE.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME, NORMAL_ALLELIC_COUNTS_FILE.getAbsolutePath())
                .addOutput(new File("Non-existent-path"))
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, outputPrefix);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = UserException.class)
    public void testNonExistentDenoisedCopyRatiosFile() {
        final File outputDir = createTempDir("testDir");
        final String outputPrefix = "test";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, "Non-existent-file")
                .addArgument(CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME, NORMAL_ALLELIC_COUNTS_FILE.getAbsolutePath())
                .addOutput(outputDir)
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, outputPrefix);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = UserException.class)
    public void testNonExistentAllelicCountsFile() {
        final File outputDir = createTempDir("testDir");
        final String outputPrefix = "test";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, TUMOR_DENOISED_COPY_RATIOS_FILE.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, "Non-existent-file")
                .addArgument(CopyNumberStandardArgument.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME, NORMAL_ALLELIC_COUNTS_FILE.getAbsolutePath())
                .addOutput(outputDir)
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, outputPrefix);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = UserException.class)
    public void testNonExistentNormalAllelicCountsFile() {
        final File outputDir = createTempDir("testDir");
        final String outputPrefix = "test";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, TUMOR_DENOISED_COPY_RATIOS_FILE.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME, "Non-existent-file")
                .addOutput(outputDir)
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, outputPrefix);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMissingAllelicCountsFile() {
        final File outputDir = createTempDir("testDir");
        final String outputPrefix = "test";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, TUMOR_DENOISED_COPY_RATIOS_FILE.getAbsolutePath())
                .addArgument(CopyNumberStandardArgument.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME, NORMAL_ALLELIC_COUNTS_FILE.getAbsolutePath())
                .addOutput(outputDir)
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, outputPrefix);
        runCommandLine(argsBuilder);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMissingCaseSampleFiles() {
        final File outputDir = createTempDir("testDir");
        final String outputPrefix = "test";
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addArgument(CopyNumberStandardArgument.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME, NORMAL_ALLELIC_COUNTS_FILE.getAbsolutePath())
                .addOutput(outputDir)
                .addArgument(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, outputPrefix);
        runCommandLine(argsBuilder);
    }
}