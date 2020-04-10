package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.*;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.models.AlleleFractionParameter;
import org.broadinstitute.hellbender.tools.copynumber.models.CopyRatioParameter;
import org.broadinstitute.hellbender.tools.copynumber.models.MultidimensionalModellerUnitTest;
import org.broadinstitute.hellbender.tools.copynumber.segmentation.MultisampleMultidimensionalKernelSegmenter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Integration tests for {@link ModelSegments}.  We test for input validation across various run modes of the tool
 * and for consistency of metadata in output, but do not test for correctness of output (which is tested elsewhere,
 * e.g. {@link MultisampleMultidimensionalKernelSegmenter} and {@link MultidimensionalModellerUnitTest}).
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class ModelSegmentsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber");
    private static final File TUMOR_1_DENOISED_COPY_RATIOS_FILE = new File(TEST_SUB_DIR,
            "model-segments-wes-tumor-1-denoised-copy-ratios-SM-74P4M-v1-chr20-downsampled.deduplicated.denoisedCR.tsv");
    private static final File TUMOR_1_DENOISED_COPY_RATIOS_WITH_MISSING_INTERVALS_FILE = new File(TEST_SUB_DIR,
            "model-segments-wes-tumor-1-denoised-copy-ratios-with-missing-intervals.denoisedCR.tsv");
    private static final File TUMOR_1_ALLELIC_COUNTS_FILE = new File(TEST_SUB_DIR,
            "model-segments-wes-tumor-1-allelic-counts-SM-74P4M-v1-chr20-downsampled.deduplicated.allelicCounts.tsv");
    private static final File TUMOR_1_ALLELIC_COUNTS_WITH_MISSING_SITES_FILE = new File(TEST_SUB_DIR,
            "model-segments-wes-tumor-1-allelic-counts-with-missing-sites.allelicCounts.tsv");
    private static final File TUMOR_2_DENOISED_COPY_RATIOS_FILE = new File(TEST_SUB_DIR,
            "model-segments-wes-tumor-2-denoised-copy-ratios-SM-74P4M-v1-chr20-downsampled.deduplicated.denoisedCR.tsv");
    private static final File TUMOR_2_ALLELIC_COUNTS_FILE = new File(TEST_SUB_DIR,
            "model-segments-wes-tumor-2-allelic-counts-SM-74P4M-v1-chr20-downsampled.deduplicated.allelicCounts.tsv");
    private static final File NORMAL_ALLELIC_COUNTS_FILE = new File(TEST_SUB_DIR,
            "model-segments-wes-normal-allelic-counts-SM-74NEG-v1-chr20-downsampled.deduplicated.allelicCounts.tsv");
    private static final File NORMAL_ALLELIC_COUNTS_WITH_MISSING_SITES_FILE = new File(TEST_SUB_DIR,
            "model-segments-wes-normal-allelic-counts-with-missing-sites.allelicCounts.tsv");
    
    private static final String OUTPUT_PREFIX = "test";
    
    private static final SampleLocatableMetadata TUMOR_1_EXPECTED_METADATA = new CopyRatioCollection(TUMOR_1_DENOISED_COPY_RATIOS_FILE).getMetadata();
    private static final SampleLocatableMetadata TUMOR_2_EXPECTED_METADATA = new CopyRatioCollection(TUMOR_2_DENOISED_COPY_RATIOS_FILE).getMetadata();
    private static final SampleLocatableMetadata NORMAL_EXPECTED_METADATA = new AllelicCountCollection(NORMAL_ALLELIC_COUNTS_FILE).getMetadata();

    @Test
    public void testMetadata() {
        Assert.assertEquals(
                TUMOR_1_EXPECTED_METADATA,
                CopyNumberArgumentValidationUtils.getValidatedMetadata(
                        new CopyRatioCollection(TUMOR_1_DENOISED_COPY_RATIOS_FILE),
                        new CopyRatioCollection(TUMOR_1_DENOISED_COPY_RATIOS_WITH_MISSING_INTERVALS_FILE),
                        new AllelicCountCollection(TUMOR_1_ALLELIC_COUNTS_FILE),
                        new AllelicCountCollection(TUMOR_1_ALLELIC_COUNTS_WITH_MISSING_SITES_FILE)));
        Assert.assertEquals(
                TUMOR_2_EXPECTED_METADATA,
                CopyNumberArgumentValidationUtils.getValidatedMetadata(
                        new CopyRatioCollection(TUMOR_2_DENOISED_COPY_RATIOS_FILE),
                        new AllelicCountCollection(TUMOR_2_ALLELIC_COUNTS_FILE)));
        Assert.assertEquals(
                NORMAL_EXPECTED_METADATA,
                CopyNumberArgumentValidationUtils.getValidatedMetadata(
                        new AllelicCountCollection(NORMAL_ALLELIC_COUNTS_FILE),
                        new AllelicCountCollection(NORMAL_ALLELIC_COUNTS_WITH_MISSING_SITES_FILE)));
    }

    @DataProvider(name = "dataValidDataModesSingleSample")
    public Object[][] dataValidDataModesSingleSample() {
        return new Object[][]{
                {
                        TUMOR_1_DENOISED_COPY_RATIOS_FILE,
                        TUMOR_1_ALLELIC_COUNTS_FILE,
                        NORMAL_ALLELIC_COUNTS_FILE
                },
                {
                        TUMOR_1_DENOISED_COPY_RATIOS_FILE,
                        TUMOR_1_ALLELIC_COUNTS_FILE,
                        null
                },
                {
                        null,
                        TUMOR_1_ALLELIC_COUNTS_FILE,
                        NORMAL_ALLELIC_COUNTS_FILE
                },
                {
                        TUMOR_1_DENOISED_COPY_RATIOS_FILE,
                        null,
                        null
                },
                {
                        null,
                        TUMOR_1_ALLELIC_COUNTS_FILE,
                        null
                }
        };
    }

    @DataProvider(name = "dataInvalidDataModesSingleSample")
    public Object[][] dataInvalidDataModesSingleSample() {
        return new Object[][]{
                //allele-fraction sites mismatch
                {
                        TUMOR_1_DENOISED_COPY_RATIOS_FILE,
                        TUMOR_1_ALLELIC_COUNTS_FILE,
                        NORMAL_ALLELIC_COUNTS_WITH_MISSING_SITES_FILE
                },
                {
                        null,
                        TUMOR_1_ALLELIC_COUNTS_FILE,
                        NORMAL_ALLELIC_COUNTS_WITH_MISSING_SITES_FILE
                },

                //sample names mismatch
                {
                        TUMOR_1_DENOISED_COPY_RATIOS_FILE,
                        TUMOR_2_ALLELIC_COUNTS_FILE,
                        NORMAL_ALLELIC_COUNTS_FILE
                },
                {
                        TUMOR_1_DENOISED_COPY_RATIOS_FILE,
                        TUMOR_2_ALLELIC_COUNTS_FILE,
                        null
                },

                //missing case allelic-counts file
                {
                        TUMOR_1_DENOISED_COPY_RATIOS_FILE,
                        null,
                        NORMAL_ALLELIC_COUNTS_FILE
                },

                //missing case files
                {
                        null,
                        null,
                        NORMAL_ALLELIC_COUNTS_FILE
                }
        };
    }

    @Test(dataProvider = "dataValidDataModesSingleSample")
    public void testValidDataModesSingleSample(final File denoisedCopyRatiosFile,
                                               final File allelicCountsFile,
                                               final File normalAllelicCountsFile) {
        final File outputDir = createTempDir("testDir");
        final ArgumentsBuilder argsBuilder = buildArgsBuilderSingleSample(
                outputDir, denoisedCopyRatiosFile, allelicCountsFile, normalAllelicCountsFile);
        runCommandLine(argsBuilder);
        final boolean isAllelicCountsPresent = allelicCountsFile != null;
        final boolean isNormalAllelicCountsPresent = normalAllelicCountsFile != null;
        assertOutputFilesSingleSample(outputDir, isAllelicCountsPresent, isNormalAllelicCountsPresent);
    }

    @Test(dataProvider = "dataInvalidDataModesSingleSample", expectedExceptions = IllegalArgumentException.class)
    public void testInvalidDataModesSingleSample(final File denoisedCopyRatiosFile,
                                                 final File allelicCountsFile,
                                                 final File normalAllelicCountsFile) {
        final File outputDir = createTempDir("testDir");
        final ArgumentsBuilder argsBuilder = buildArgsBuilderSingleSample(
                outputDir, denoisedCopyRatiosFile, allelicCountsFile, normalAllelicCountsFile);
        runCommandLine(argsBuilder);
    }

    private static ArgumentsBuilder buildArgsBuilderSingleSample(final File outputDir,
                                                                 final File denoisedCopyRatiosFile,
                                                                 final File allelicCountsFile,
                                                                 final File normalAllelicCountsFile) {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addOutput(outputDir)
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX);
        if (denoisedCopyRatiosFile != null) {
            argsBuilder.add(CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, denoisedCopyRatiosFile);
        }
        if (allelicCountsFile != null) {
            argsBuilder.add(CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, allelicCountsFile);
        }
        if (normalAllelicCountsFile != null) {
            argsBuilder.add(CopyNumberStandardArgument.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME, normalAllelicCountsFile);
        }
        return argsBuilder;
    }

    private static void assertOutputFilesSingleSample(final File outputDir,
                                                      final boolean isAllelicCountsPresent,
                                                      final boolean isNormalAllelicCountsPresent) {
        Assert.assertFalse(!isAllelicCountsPresent && isNormalAllelicCountsPresent);
        for (final String fileTag : Arrays.asList(ModelSegments.BEGIN_FIT_FILE_TAG, ModelSegments.FINAL_FIT_FILE_TAG)) {
            final ModeledSegmentCollection modeledSegments = new ModeledSegmentCollection(
                    new File(outputDir, OUTPUT_PREFIX + fileTag + ModelSegments.SEGMENTS_FILE_SUFFIX));
            Assert.assertEquals(TUMOR_1_EXPECTED_METADATA, modeledSegments.getMetadata());

            final ParameterDecileCollection<CopyRatioParameter> copyRatioParameters = new ParameterDecileCollection<>(
                    new File(outputDir, OUTPUT_PREFIX + fileTag + ModelSegments.COPY_RATIO_MODEL_PARAMETER_FILE_SUFFIX), CopyRatioParameter.class);
            Assert.assertEquals(TUMOR_1_EXPECTED_METADATA.getSampleName(), copyRatioParameters.getMetadata().getSampleName());

            final ParameterDecileCollection<AlleleFractionParameter> alleleFractionParameters = new ParameterDecileCollection<>(
                    new File(outputDir, OUTPUT_PREFIX + fileTag + ModelSegments.ALLELE_FRACTION_MODEL_PARAMETER_FILE_SUFFIX), AlleleFractionParameter.class);
            Assert.assertEquals(TUMOR_1_EXPECTED_METADATA.getSampleName(), alleleFractionParameters.getMetadata().getSampleName());
        }

        final CopyRatioSegmentCollection copyRatioSegments = new CopyRatioSegmentCollection(
                new File(outputDir, OUTPUT_PREFIX + ModelSegments.COPY_RATIO_SEGMENTS_FOR_CALLER_FILE_SUFFIX));
        Assert.assertEquals(TUMOR_1_EXPECTED_METADATA, copyRatioSegments.getMetadata());

        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + ModelSegments.COPY_RATIO_LEGACY_SEGMENTS_FILE_SUFFIX).exists());
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + ModelSegments.ALLELE_FRACTION_LEGACY_SEGMENTS_FILE_SUFFIX).exists());

        AllelicCountCollection hetAllelicCounts = null;
        if (isAllelicCountsPresent) {
            hetAllelicCounts = new AllelicCountCollection(
                    new File(outputDir, OUTPUT_PREFIX + ModelSegments.HET_ALLELIC_COUNTS_FILE_SUFFIX));
            Assert.assertEquals(TUMOR_1_EXPECTED_METADATA, hetAllelicCounts.getMetadata());
        }
        if (isNormalAllelicCountsPresent) { //if this is true, case sample allelic counts will be present
            final AllelicCountCollection hetNormalAllelicCounts = new AllelicCountCollection(
                    new File(outputDir, OUTPUT_PREFIX + ModelSegments.NORMAL_HET_ALLELIC_COUNTS_FILE_SUFFIX));
            Assert.assertEquals(NORMAL_EXPECTED_METADATA, hetNormalAllelicCounts.getMetadata());
            Assert.assertTrue(CopyNumberArgumentValidationUtils.isSameDictionary(               //sequence dictionary should be the same
                    TUMOR_1_EXPECTED_METADATA.getSequenceDictionary(), hetNormalAllelicCounts.getMetadata().getSequenceDictionary()));
            Assert.assertEquals(hetAllelicCounts.getIntervals(), hetNormalAllelicCounts.getIntervals());
        }

        Assert.assertFalse(new File(outputDir, OUTPUT_PREFIX + ModelSegments.PICARD_INTERVAL_LIST_FILE_SUFFIX).exists());
    }

    @DataProvider(name = "dataValidDataModesMultipleSamples")
    public Object[][] dataValidDataModesMultipleSamples() {
        return new Object[][]{
                {
                        Arrays.asList(TUMOR_1_DENOISED_COPY_RATIOS_FILE, TUMOR_2_DENOISED_COPY_RATIOS_FILE),
                        Arrays.asList(TUMOR_1_ALLELIC_COUNTS_FILE, TUMOR_2_ALLELIC_COUNTS_FILE),
                        NORMAL_ALLELIC_COUNTS_FILE
                },
                {
                        Arrays.asList(TUMOR_1_DENOISED_COPY_RATIOS_FILE, TUMOR_2_DENOISED_COPY_RATIOS_FILE),
                        Arrays.asList(TUMOR_1_ALLELIC_COUNTS_FILE, TUMOR_2_ALLELIC_COUNTS_FILE),
                        null
                },
                {
                        null,
                        Arrays.asList(TUMOR_1_ALLELIC_COUNTS_FILE, TUMOR_2_ALLELIC_COUNTS_FILE),
                        NORMAL_ALLELIC_COUNTS_FILE
                },
                {
                        Arrays.asList(TUMOR_1_DENOISED_COPY_RATIOS_FILE, TUMOR_2_DENOISED_COPY_RATIOS_FILE),
                        null,
                        null
                },
                {
                        null,
                        Arrays.asList(TUMOR_1_ALLELIC_COUNTS_FILE, TUMOR_2_ALLELIC_COUNTS_FILE),
                        null
                }
        };
    }

    @DataProvider(name = "dataInvalidDataModesMultipleSamples")
    public Object[][] dataInvalidDataModesMultipleSamples() {
        return new Object[][]{
                //copy-ratio intervals mismatch
                {
                        Arrays.asList(TUMOR_1_DENOISED_COPY_RATIOS_WITH_MISSING_INTERVALS_FILE, TUMOR_2_DENOISED_COPY_RATIOS_FILE),
                        Arrays.asList(TUMOR_1_ALLELIC_COUNTS_FILE, TUMOR_2_ALLELIC_COUNTS_FILE),
                        NORMAL_ALLELIC_COUNTS_FILE
                },
                {
                        Arrays.asList(TUMOR_1_DENOISED_COPY_RATIOS_WITH_MISSING_INTERVALS_FILE, TUMOR_2_DENOISED_COPY_RATIOS_FILE),
                        Arrays.asList(TUMOR_1_ALLELIC_COUNTS_FILE, TUMOR_2_ALLELIC_COUNTS_FILE),
                        null
                },
                {
                        Arrays.asList(TUMOR_1_DENOISED_COPY_RATIOS_WITH_MISSING_INTERVALS_FILE, TUMOR_2_DENOISED_COPY_RATIOS_FILE),
                        null,
                        null
                },

                //case allele-fraction sites mismatch
                {
                        Arrays.asList(TUMOR_1_DENOISED_COPY_RATIOS_FILE, TUMOR_2_DENOISED_COPY_RATIOS_FILE),
                        Arrays.asList(TUMOR_1_ALLELIC_COUNTS_WITH_MISSING_SITES_FILE, TUMOR_2_ALLELIC_COUNTS_FILE),
                        NORMAL_ALLELIC_COUNTS_FILE
                },
                {
                        Arrays.asList(TUMOR_1_DENOISED_COPY_RATIOS_FILE, TUMOR_2_DENOISED_COPY_RATIOS_FILE),
                        Arrays.asList(TUMOR_1_ALLELIC_COUNTS_WITH_MISSING_SITES_FILE, TUMOR_2_ALLELIC_COUNTS_FILE),
                        null
                },
                {
                        null,
                        Arrays.asList(TUMOR_1_ALLELIC_COUNTS_WITH_MISSING_SITES_FILE, TUMOR_2_ALLELIC_COUNTS_FILE),
                        null
                },

                //normal allele-fraction sites mismatch
                {
                        Arrays.asList(TUMOR_1_DENOISED_COPY_RATIOS_FILE, TUMOR_2_DENOISED_COPY_RATIOS_FILE),
                        Arrays.asList(TUMOR_1_ALLELIC_COUNTS_FILE, TUMOR_2_ALLELIC_COUNTS_FILE),
                        NORMAL_ALLELIC_COUNTS_WITH_MISSING_SITES_FILE
                },
                {
                        null,
                        Arrays.asList(TUMOR_1_ALLELIC_COUNTS_FILE, TUMOR_2_ALLELIC_COUNTS_FILE),
                        NORMAL_ALLELIC_COUNTS_WITH_MISSING_SITES_FILE
                },

                //sample order mismatch
                {
                        Arrays.asList(TUMOR_1_DENOISED_COPY_RATIOS_FILE, TUMOR_2_DENOISED_COPY_RATIOS_FILE),
                        Arrays.asList(TUMOR_2_ALLELIC_COUNTS_FILE, TUMOR_1_ALLELIC_COUNTS_FILE),
                        NORMAL_ALLELIC_COUNTS_FILE
                },
                {
                        Arrays.asList(TUMOR_1_DENOISED_COPY_RATIOS_FILE, TUMOR_2_DENOISED_COPY_RATIOS_FILE),
                        Arrays.asList(TUMOR_2_ALLELIC_COUNTS_FILE, TUMOR_1_ALLELIC_COUNTS_FILE),
                        null
                },

                //sample number mismatch
                {
                        Arrays.asList(TUMOR_1_DENOISED_COPY_RATIOS_FILE, TUMOR_2_DENOISED_COPY_RATIOS_FILE),
                        Collections.singletonList(TUMOR_1_ALLELIC_COUNTS_FILE),
                        NORMAL_ALLELIC_COUNTS_FILE
                },
                {
                        Arrays.asList(TUMOR_1_DENOISED_COPY_RATIOS_FILE, TUMOR_2_DENOISED_COPY_RATIOS_FILE),
                        Collections.singletonList(TUMOR_1_ALLELIC_COUNTS_FILE),
                        null
                },
                {
                        Collections.singletonList(TUMOR_1_DENOISED_COPY_RATIOS_FILE),
                        Arrays.asList(TUMOR_1_ALLELIC_COUNTS_FILE, TUMOR_2_ALLELIC_COUNTS_FILE),
                        NORMAL_ALLELIC_COUNTS_FILE
                },
                {
                        Collections.singletonList(TUMOR_1_DENOISED_COPY_RATIOS_FILE),
                        Arrays.asList(TUMOR_1_ALLELIC_COUNTS_FILE, TUMOR_2_ALLELIC_COUNTS_FILE),
                        null
                },

                //missing case allelic-counts files
                {
                        Arrays.asList(TUMOR_1_DENOISED_COPY_RATIOS_FILE, TUMOR_2_DENOISED_COPY_RATIOS_FILE),
                        null,
                        NORMAL_ALLELIC_COUNTS_FILE
                }
        };
    }

    @Test(dataProvider = "dataValidDataModesMultipleSamples")
    public void testValidDataModesMultipleSamples(final List<File> denoisedCopyRatiosFiles,
                                                  final List<File> allelicCountsFiles,
                                                  final File normalAllelicCountsFile) {
        final File outputDir = createTempDir("testDir");
        final ArgumentsBuilder argsBuilder = buildArgsBuilderMultipleSamples(
                outputDir, denoisedCopyRatiosFiles, allelicCountsFiles, normalAllelicCountsFile);
        runCommandLine(argsBuilder);
        final boolean isAllelicCountsPresent = allelicCountsFiles != null;
        final boolean isNormalAllelicCountsPresent = normalAllelicCountsFile != null;
        assertOutputFilesMultipleSamples(outputDir, isAllelicCountsPresent, isNormalAllelicCountsPresent);
    }

    @Test(dataProvider = "dataInvalidDataModesMultipleSamples", expectedExceptions = IllegalArgumentException.class)
    public void testInvalidDataModesMultipleSamples(final List<File> denoisedCopyRatiosFiles,
                                                    final List<File> allelicCountsFiles,
                                                    final File normalAllelicCountsFile) {
        final File outputDir = createTempDir("testDir");
        final ArgumentsBuilder argsBuilder = buildArgsBuilderMultipleSamples(
                outputDir, denoisedCopyRatiosFiles, allelicCountsFiles, normalAllelicCountsFile);
        runCommandLine(argsBuilder);
    }

    private static ArgumentsBuilder buildArgsBuilderMultipleSamples(final File outputDir,
                                                                    final List<File> denoisedCopyRatiosFiles,
                                                                    final List<File> allelicCountsFiles,
                                                                    final File normalAllelicCountsFile) {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addOutput(outputDir)
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX);
        if (denoisedCopyRatiosFiles != null) {
            denoisedCopyRatiosFiles.forEach(f -> argsBuilder.add(CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, f));
        }
        if (allelicCountsFiles != null) {
            allelicCountsFiles.forEach(f -> argsBuilder.add(CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, f));
        }
        if (normalAllelicCountsFile != null) {
            argsBuilder.add(CopyNumberStandardArgument.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME, normalAllelicCountsFile);
        }
        return argsBuilder;
    }

    private static void assertOutputFilesMultipleSamples(final File outputDir,
                                                         final boolean isAllelicCountsPresent,
                                                         final boolean isNormalAllelicCountsPresent) {
        Assert.assertFalse(!isAllelicCountsPresent && isNormalAllelicCountsPresent);
        for (final String fileTag : Arrays.asList(ModelSegments.BEGIN_FIT_FILE_TAG, ModelSegments.FINAL_FIT_FILE_TAG)) {
            Assert.assertFalse(new File(outputDir, OUTPUT_PREFIX + fileTag + ModelSegments.SEGMENTS_FILE_SUFFIX).exists());
            Assert.assertFalse(new File(outputDir, OUTPUT_PREFIX + fileTag + ModelSegments.COPY_RATIO_MODEL_PARAMETER_FILE_SUFFIX).exists());
            Assert.assertFalse(new File(outputDir, OUTPUT_PREFIX + fileTag + ModelSegments.ALLELE_FRACTION_MODEL_PARAMETER_FILE_SUFFIX).exists());
        }
        Assert.assertFalse(new File(outputDir, OUTPUT_PREFIX + ModelSegments.COPY_RATIO_SEGMENTS_FOR_CALLER_FILE_SUFFIX).exists());
        Assert.assertFalse(new File(outputDir, OUTPUT_PREFIX + ModelSegments.COPY_RATIO_LEGACY_SEGMENTS_FILE_SUFFIX).exists());
        Assert.assertFalse(new File(outputDir, OUTPUT_PREFIX + ModelSegments.ALLELE_FRACTION_LEGACY_SEGMENTS_FILE_SUFFIX).exists());
        Assert.assertFalse(new File(outputDir, OUTPUT_PREFIX + ModelSegments.HET_ALLELIC_COUNTS_FILE_SUFFIX).exists());
        Assert.assertFalse(new File(outputDir, OUTPUT_PREFIX + ModelSegments.NORMAL_HET_ALLELIC_COUNTS_FILE_SUFFIX).exists());

        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + ModelSegments.PICARD_INTERVAL_LIST_FILE_SUFFIX).exists());
    }
}
