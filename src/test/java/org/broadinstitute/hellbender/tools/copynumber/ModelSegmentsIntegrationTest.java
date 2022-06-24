package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.CopyNumberTestUtils;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.ModeledSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.ParameterDecileCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.models.AlleleFractionParameter;
import org.broadinstitute.hellbender.tools.copynumber.models.CopyRatioParameter;
import org.broadinstitute.hellbender.tools.copynumber.models.MultidimensionalModellerUnitTest;
import org.broadinstitute.hellbender.tools.copynumber.segmentation.MultisampleMultidimensionalKernelSegmenter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Integration tests for {@link ModelSegments}.  We test for input validation across various run modes of the tool
 * and for consistency of metadata in output, but do not test for correctness of output (which is tested elsewhere,
 * e.g. {@link MultisampleMultidimensionalKernelSegmenter} and {@link MultidimensionalModellerUnitTest}).
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class ModelSegmentsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber");
    private static final File LARGE_TEST_RESOURCES_SUB_DIR = new File(largeFileTestDir, "org/broadinstitute/hellbender/tools/copynumber");
    private static final File EXACT_MATCH_EXPECTED_SUB_DIR = new File(TEST_SUB_DIR, "model-segments-expected");
    private static final File TUMOR_1_DENOISED_COPY_RATIOS_FILE = new File(LARGE_TEST_RESOURCES_SUB_DIR,
            "model-segments-wes-tumor-1-denoised-copy-ratios-SM-74P4M-v1-chr20.denoisedCR.tsv");
    private static final File TUMOR_1_DENOISED_COPY_RATIOS_WITH_MISSING_INTERVALS_FILE = new File(LARGE_TEST_RESOURCES_SUB_DIR,
            "model-segments-wes-tumor-1-denoised-copy-ratios-with-missing-intervals.denoisedCR.tsv");
    private static final File TUMOR_1_ALLELIC_COUNTS_FILE = new File(LARGE_TEST_RESOURCES_SUB_DIR,
            "model-segments-wes-tumor-1-allelic-counts-SM-74P4M-v1-chr20.allelicCounts.tsv");
    private static final File TUMOR_1_ALLELIC_COUNTS_WITH_MISSING_SITES_FILE = new File(LARGE_TEST_RESOURCES_SUB_DIR,
            "model-segments-wes-tumor-1-allelic-counts-with-missing-sites.allelicCounts.tsv");
    private static final File TUMOR_2_DENOISED_COPY_RATIOS_FILE = new File(LARGE_TEST_RESOURCES_SUB_DIR,
            "model-segments-wes-tumor-2-denoised-copy-ratios-SM-74P4M-v1-chr20.denoisedCR.tsv");
    private static final File TUMOR_2_ALLELIC_COUNTS_FILE = new File(LARGE_TEST_RESOURCES_SUB_DIR,
            "model-segments-wes-tumor-2-allelic-counts-SM-74P4M-v1-chr20.allelicCounts.tsv");
    private static final File NORMAL_ALLELIC_COUNTS_FILE = new File(LARGE_TEST_RESOURCES_SUB_DIR,
            "model-segments-wes-normal-allelic-counts-SM-74NEG-v1-chr20.allelicCounts.tsv");
    private static final File NORMAL_ALLELIC_COUNTS_WITH_MISSING_SITES_FILE = new File(LARGE_TEST_RESOURCES_SUB_DIR,
            "model-segments-wes-normal-allelic-counts-with-missing-sites.allelicCounts.tsv");
    
    private static final String DEFAULT_OUTPUT_PREFIX = "test";
    
    private static final SampleLocatableMetadata TUMOR_1_EXPECTED_METADATA = new CopyRatioCollection(TUMOR_1_DENOISED_COPY_RATIOS_FILE).getMetadata();
    private static final SampleLocatableMetadata TUMOR_2_EXPECTED_METADATA = new CopyRatioCollection(TUMOR_2_DENOISED_COPY_RATIOS_FILE).getMetadata();
    private static final SampleLocatableMetadata NORMAL_EXPECTED_METADATA = new AllelicCountCollection(NORMAL_ALLELIC_COUNTS_FILE).getMetadata();

    /**
     * Note that {@link org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberFormatsUtils#DOUBLE_FORMAT}
     * is set so that doubles in somatic CNV outputs will have 6 decimal places. We thus set the allowed delta
     * to detect differences at that level.
     */
    private static final double ALLOWED_DELTA_FOR_DOUBLE_VALUES = 1E-6;

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
                    "single-sample-cr-ac-nac", // cr-ac-nac = copy ratios + allelic counts + normal allelic counts
                    TUMOR_1_DENOISED_COPY_RATIOS_FILE,
                    TUMOR_1_ALLELIC_COUNTS_FILE,
                    NORMAL_ALLELIC_COUNTS_FILE
                },
                {
                    "single-sample-cr-ac",
                    TUMOR_1_DENOISED_COPY_RATIOS_FILE,
                    TUMOR_1_ALLELIC_COUNTS_FILE,
                    null
                },
                {
                    "single-sample-ac-nac",
                    null,
                    TUMOR_1_ALLELIC_COUNTS_FILE,
                    NORMAL_ALLELIC_COUNTS_FILE
                },
                {
                    "single-sample-cr",
                    TUMOR_1_DENOISED_COPY_RATIOS_FILE,
                    null,
                    null
                },
                {
                    "single-sample-ac",
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
    public void testValidDataModesSingleSample(final String outputPrefix,
                                               final File denoisedCopyRatiosFile,
                                               final File allelicCountsFile,
                                               final File normalAllelicCountsFile) {
        final File outputDir = createTempDir("testDir");
        final ArgumentsBuilder argsBuilder = buildArgsBuilderSingleSample(
                outputDir, outputPrefix, denoisedCopyRatiosFile, allelicCountsFile, normalAllelicCountsFile);
        runCommandLine(argsBuilder);
        final boolean isAllelicCountsPresent = allelicCountsFile != null;
        final boolean isNormalAllelicCountsPresent = normalAllelicCountsFile != null;
        assertOutputFilesSingleSample(outputDir, outputPrefix, isAllelicCountsPresent, isNormalAllelicCountsPresent);
    }

    @Test(dataProvider = "dataInvalidDataModesSingleSample", expectedExceptions = IllegalArgumentException.class)
    public void testInvalidDataModesSingleSample(final File denoisedCopyRatiosFile,
                                                 final File allelicCountsFile,
                                                 final File normalAllelicCountsFile) {
        final File outputDir = createTempDir("testDir");
        final ArgumentsBuilder argsBuilder = buildArgsBuilderSingleSample(
                outputDir, DEFAULT_OUTPUT_PREFIX, denoisedCopyRatiosFile, allelicCountsFile, normalAllelicCountsFile);
        runCommandLine(argsBuilder);
    }

    private static ArgumentsBuilder buildArgsBuilderSingleSample(final File outputDir,
                                                                 final String outputPrefix,
                                                                 final File denoisedCopyRatiosFile,
                                                                 final File allelicCountsFile,
                                                                 final File normalAllelicCountsFile) {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addOutput(outputDir)
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, outputPrefix);
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

    /**
     * Note that exact-match tests do not completely cover reading/parsing of output files;
     * we assume that it is unlikely that reading/parsing code will change/break and that rough coverage in e.g.
     * {@link org.broadinstitute.hellbender.tools.copynumber.formats.collections.AbstractSampleLocatableCollectionUnitTest}
     * will suffice.
     */
    private static void assertOutputFilesSingleSample(final File outputDir,
                                                      final String outputPrefix,
                                                      final boolean isAllelicCountsPresent,
                                                      final boolean isNormalAllelicCountsPresent) {
        Assert.assertFalse(!isAllelicCountsPresent && isNormalAllelicCountsPresent);
        for (final String fileTag : Arrays.asList(ModelSegments.BEGIN_FIT_FILE_TAG, ModelSegments.FINAL_FIT_FILE_TAG)) {
            // test exact match with expected outputs
            for (final String suffix : Arrays.asList(
                    ModelSegments.SEGMENTS_FILE_SUFFIX,
                    ModelSegments.COPY_RATIO_MODEL_PARAMETER_FILE_SUFFIX,
                    ModelSegments.ALLELE_FRACTION_MODEL_PARAMETER_FILE_SUFFIX)) {
                CopyNumberTestUtils.assertFilesEqualUpToAllowedDeltaForDoubleValues(
                        new File(outputDir, outputPrefix + fileTag + suffix),
                        new File(EXACT_MATCH_EXPECTED_SUB_DIR, outputPrefix + fileTag + suffix),
                        ALLOWED_DELTA_FOR_DOUBLE_VALUES,
                        logger);
            }

            // test sample-name reading
            final ModeledSegmentCollection modeledSegments = new ModeledSegmentCollection(
                    new File(outputDir, outputPrefix + fileTag + ModelSegments.SEGMENTS_FILE_SUFFIX));
            Assert.assertEquals(TUMOR_1_EXPECTED_METADATA, modeledSegments.getMetadata());

            final ParameterDecileCollection<CopyRatioParameter> copyRatioParameters = new ParameterDecileCollection<>(
                    new File(outputDir, outputPrefix + fileTag + ModelSegments.COPY_RATIO_MODEL_PARAMETER_FILE_SUFFIX), CopyRatioParameter.class);
            Assert.assertEquals(TUMOR_1_EXPECTED_METADATA.getSampleName(), copyRatioParameters.getMetadata().getSampleName());

            final ParameterDecileCollection<AlleleFractionParameter> alleleFractionParameters = new ParameterDecileCollection<>(
                    new File(outputDir, outputPrefix + fileTag + ModelSegments.ALLELE_FRACTION_MODEL_PARAMETER_FILE_SUFFIX), AlleleFractionParameter.class);
            Assert.assertEquals(TUMOR_1_EXPECTED_METADATA.getSampleName(), alleleFractionParameters.getMetadata().getSampleName());
        }

        CopyNumberTestUtils.assertFilesEqualUpToAllowedDeltaForDoubleValues(
                new File(outputDir, outputPrefix  + ModelSegments.COPY_RATIO_SEGMENTS_FOR_CALLER_FILE_SUFFIX),
                new File(EXACT_MATCH_EXPECTED_SUB_DIR, outputPrefix  + ModelSegments.COPY_RATIO_SEGMENTS_FOR_CALLER_FILE_SUFFIX),
                ALLOWED_DELTA_FOR_DOUBLE_VALUES,
                logger);
        final CopyRatioSegmentCollection copyRatioSegments = new CopyRatioSegmentCollection(
                new File(outputDir, outputPrefix + ModelSegments.COPY_RATIO_SEGMENTS_FOR_CALLER_FILE_SUFFIX));
        Assert.assertEquals(TUMOR_1_EXPECTED_METADATA, copyRatioSegments.getMetadata());

        CopyNumberTestUtils.assertFilesEqualUpToAllowedDeltaForDoubleValues(
                new File(outputDir, outputPrefix  + ModelSegments.COPY_RATIO_LEGACY_SEGMENTS_FILE_SUFFIX),
                new File(EXACT_MATCH_EXPECTED_SUB_DIR, outputPrefix  + ModelSegments.COPY_RATIO_LEGACY_SEGMENTS_FILE_SUFFIX),
                ALLOWED_DELTA_FOR_DOUBLE_VALUES,
                logger);
        CopyNumberTestUtils.assertFilesEqualUpToAllowedDeltaForDoubleValues(
                new File(outputDir, outputPrefix  + ModelSegments.ALLELE_FRACTION_LEGACY_SEGMENTS_FILE_SUFFIX),
                new File(EXACT_MATCH_EXPECTED_SUB_DIR, outputPrefix  + ModelSegments.ALLELE_FRACTION_LEGACY_SEGMENTS_FILE_SUFFIX),
                ALLOWED_DELTA_FOR_DOUBLE_VALUES,
                logger);

        AllelicCountCollection hetAllelicCounts = null;
        if (isAllelicCountsPresent) {
            CopyNumberTestUtils.assertFilesEqualUpToAllowedDeltaForDoubleValues(
                    new File(outputDir, outputPrefix  + ModelSegments.HET_ALLELIC_COUNTS_FILE_SUFFIX),
                    new File(EXACT_MATCH_EXPECTED_SUB_DIR, outputPrefix  + ModelSegments.HET_ALLELIC_COUNTS_FILE_SUFFIX),
                    ALLOWED_DELTA_FOR_DOUBLE_VALUES,
                    logger);
            hetAllelicCounts = new AllelicCountCollection(
                    new File(outputDir, outputPrefix + ModelSegments.HET_ALLELIC_COUNTS_FILE_SUFFIX));
            Assert.assertEquals(TUMOR_1_EXPECTED_METADATA, hetAllelicCounts.getMetadata());
        }
        if (isNormalAllelicCountsPresent) { //if this is true, case sample allelic counts will be present
            CopyNumberTestUtils.assertFilesEqualUpToAllowedDeltaForDoubleValues(
                    new File(outputDir, outputPrefix  + ModelSegments.NORMAL_HET_ALLELIC_COUNTS_FILE_SUFFIX),
                    new File(EXACT_MATCH_EXPECTED_SUB_DIR, outputPrefix  + ModelSegments.NORMAL_HET_ALLELIC_COUNTS_FILE_SUFFIX),
                    ALLOWED_DELTA_FOR_DOUBLE_VALUES,
                    logger);
            final AllelicCountCollection hetNormalAllelicCounts = new AllelicCountCollection(
                    new File(outputDir, outputPrefix + ModelSegments.NORMAL_HET_ALLELIC_COUNTS_FILE_SUFFIX));
            Assert.assertEquals(NORMAL_EXPECTED_METADATA, hetNormalAllelicCounts.getMetadata());
            Assert.assertTrue(CopyNumberArgumentValidationUtils.isSameDictionary(               //sequence dictionary should be the same
                    TUMOR_1_EXPECTED_METADATA.getSequenceDictionary(), hetNormalAllelicCounts.getMetadata().getSequenceDictionary()));
            Assert.assertEquals(hetAllelicCounts.getIntervals(), hetNormalAllelicCounts.getIntervals());
        }

        Assert.assertFalse(new File(outputDir, outputPrefix + ModelSegments.PICARD_INTERVAL_LIST_FILE_SUFFIX).exists());
    }

    @DataProvider(name = "dataValidDataModesMultipleSamples")
    public Object[][] dataValidDataModesMultipleSamples() {
        return new Object[][]{
                {
                    "multiple-sample-cr-ac-nac", // cr-ac-nac = copy ratios + allelic counts + normal allelic counts
                    Arrays.asList(TUMOR_1_DENOISED_COPY_RATIOS_FILE, TUMOR_2_DENOISED_COPY_RATIOS_FILE),
                    Arrays.asList(TUMOR_1_ALLELIC_COUNTS_FILE, TUMOR_2_ALLELIC_COUNTS_FILE),
                    NORMAL_ALLELIC_COUNTS_FILE
                },
                {
                    "multiple-sample-cr-ac",
                    Arrays.asList(TUMOR_1_DENOISED_COPY_RATIOS_FILE, TUMOR_2_DENOISED_COPY_RATIOS_FILE),
                    Arrays.asList(TUMOR_1_ALLELIC_COUNTS_FILE, TUMOR_2_ALLELIC_COUNTS_FILE),
                    null
                },
                {
                    "multiple-sample-ac-nac",
                    null,
                    Arrays.asList(TUMOR_1_ALLELIC_COUNTS_FILE, TUMOR_2_ALLELIC_COUNTS_FILE),
                    NORMAL_ALLELIC_COUNTS_FILE
                },
                {
                    "multiple-sample-cr",
                    Arrays.asList(TUMOR_1_DENOISED_COPY_RATIOS_FILE, TUMOR_2_DENOISED_COPY_RATIOS_FILE),
                    null,
                    null
                },
                {
                    "multiple-sample-ac",
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
    public void testValidDataModesMultipleSamples(final String outputPrefix,
                                                  final List<File> denoisedCopyRatiosFiles,
                                                  final List<File> allelicCountsFiles,
                                                  final File normalAllelicCountsFile) {
        final File outputDir = createTempDir("testDir");

        // test joint segmentation
        final ArgumentsBuilder argsBuilder = buildArgsBuilderMultipleSamples(
                outputDir, outputPrefix, denoisedCopyRatiosFiles, allelicCountsFiles, normalAllelicCountsFile);
        runCommandLine(argsBuilder);
        final boolean isAllelicCountsPresent = allelicCountsFiles != null;
        final boolean isNormalAllelicCountsPresent = normalAllelicCountsFile != null;
        assertOutputFilesMultipleSamples(outputDir, outputPrefix, isAllelicCountsPresent, isNormalAllelicCountsPresent);

        // test using joint segmentation as input to scatter of first case sample
        final ArgumentsBuilder argsBuilderSingleSample = buildArgsBuilderSingleSample(
                outputDir, outputPrefix + "-tumor-1",
                denoisedCopyRatiosFiles == null ? null : denoisedCopyRatiosFiles.get(0),
                allelicCountsFiles == null ? null : allelicCountsFiles.get(0),
                normalAllelicCountsFile);
        argsBuilderSingleSample.add(    // add the joint segmentation as input
                CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME,
                new File(outputDir, outputPrefix + ModelSegments.PICARD_INTERVAL_LIST_FILE_SUFFIX));
        runCommandLine(argsBuilderSingleSample);
        assertOutputFilesSingleSample(outputDir, outputPrefix + "-tumor-1", isAllelicCountsPresent, isNormalAllelicCountsPresent);
    }

    @Test(dataProvider = "dataInvalidDataModesMultipleSamples", expectedExceptions = IllegalArgumentException.class)
    public void testInvalidDataModesMultipleSamples(final List<File> denoisedCopyRatiosFiles,
                                                    final List<File> allelicCountsFiles,
                                                    final File normalAllelicCountsFile) {
        final File outputDir = createTempDir("testDir");

        // test joint segmentation
        final ArgumentsBuilder argsBuilder = buildArgsBuilderMultipleSamples(
                outputDir, DEFAULT_OUTPUT_PREFIX, denoisedCopyRatiosFiles, allelicCountsFiles, normalAllelicCountsFile);
        runCommandLine(argsBuilder);
    }

    private static ArgumentsBuilder buildArgsBuilderMultipleSamples(final File outputDir,
                                                                    final String outputPrefix,
                                                                    final List<File> denoisedCopyRatiosFiles,
                                                                    final List<File> allelicCountsFiles,
                                                                    final File normalAllelicCountsFile) {
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addOutput(outputDir)
                .add(CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, outputPrefix);
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
                                                         final String outputPrefix,
                                                         final boolean isAllelicCountsPresent,
                                                         final boolean isNormalAllelicCountsPresent) {
        Assert.assertFalse(!isAllelicCountsPresent && isNormalAllelicCountsPresent);
        for (final String fileTag : Arrays.asList(ModelSegments.BEGIN_FIT_FILE_TAG, ModelSegments.FINAL_FIT_FILE_TAG)) {
            for (final String suffix : Arrays.asList(
                    ModelSegments.SEGMENTS_FILE_SUFFIX,
                    ModelSegments.COPY_RATIO_MODEL_PARAMETER_FILE_SUFFIX,
                    ModelSegments.ALLELE_FRACTION_MODEL_PARAMETER_FILE_SUFFIX)) {
                Assert.assertFalse(new File(outputDir, outputPrefix + fileTag + suffix).exists());
            }
        }
        Assert.assertFalse(new File(outputDir, outputPrefix + ModelSegments.COPY_RATIO_SEGMENTS_FOR_CALLER_FILE_SUFFIX).exists());
        Assert.assertFalse(new File(outputDir, outputPrefix + ModelSegments.COPY_RATIO_LEGACY_SEGMENTS_FILE_SUFFIX).exists());
        Assert.assertFalse(new File(outputDir, outputPrefix + ModelSegments.ALLELE_FRACTION_LEGACY_SEGMENTS_FILE_SUFFIX).exists());
        Assert.assertFalse(new File(outputDir, outputPrefix + ModelSegments.HET_ALLELIC_COUNTS_FILE_SUFFIX).exists());
        Assert.assertFalse(new File(outputDir, outputPrefix + ModelSegments.NORMAL_HET_ALLELIC_COUNTS_FILE_SUFFIX).exists());

        CopyNumberTestUtils.assertFilesEqualUpToAllowedDeltaForDoubleValues(
                new File(outputDir, outputPrefix + ModelSegments.PICARD_INTERVAL_LIST_FILE_SUFFIX),
                new File(EXACT_MATCH_EXPECTED_SUB_DIR, outputPrefix + ModelSegments.PICARD_INTERVAL_LIST_FILE_SUFFIX),
                ALLOWED_DELTA_FOR_DOUBLE_VALUES,
                logger);
    }
}
