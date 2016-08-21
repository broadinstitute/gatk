package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.samplenamefinder.SampleNameFinder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

import static org.broadinstitute.hellbender.tools.exome.AllelicCNV.*;

/**
 * Integration tests for {@link AllelicCNV}.  Note that behavior of the tests depends on the content of the input files
 * and that changing them may break the tests.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class AllelicCNVIntegrationTest extends CommandLineProgramTest {
    private static final String DUMMY_OUTPUT_PREFIX = "dummyOutputPrefix";
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome";
    private static final File COVERAGES_FILE = new File(TEST_SUB_DIR, "coverages-for-allelic-integration.tsv");
    private static final File TUMOR_ALLELIC_COUNTS_FILE = new File(TEST_SUB_DIR, "snps-for-allelic-integration.tsv");
    private static final File SEGMENT_FILE = new File(TEST_SUB_DIR, "segments-for-allelic-integration.seg");
    private static final File ALLELIC_PON_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-pon-normal.tsv");
    private static final String SAMPLE_NAME = "test";

    @DataProvider(name = "dataACNVArgumentValidation")
    public Object[][] dataACNVArgumentValidation() {
        return new Object[][]{
                {new String[]{"--" + MAX_NUM_SNP_SEGMENTATION_ITERATIONS_LONG_NAME, "0"}},
                {new String[]{"--" + MAX_NUM_SNP_SEGMENTATION_ITERATIONS_LONG_NAME, "-1"}},
                {new String[]{"--" + SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD_LONG_NAME, "-1"}},
                {new String[]{"--" + NUM_SAMPLES_COPY_RATIO_LONG_NAME, "0"}},
                {new String[]{"--" + NUM_SAMPLES_COPY_RATIO_LONG_NAME, "-1"}},
                {new String[]{"--" + NUM_SAMPLES_COPY_RATIO_LONG_NAME, "5",
                              "--" + NUM_BURN_IN_COPY_RATIO_LONG_NAME, "10"}},
                {new String[]{"--" + NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME, "0"}},
                {new String[]{"--" + NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME, "-1"}},
                {new String[]{"--" + NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME, "5",
                              "--" + NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME, "10"}},
                {new String[]{"--" + INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME, "0"}},
                {new String[]{"--" + INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME, "-1"}},
                {new String[]{"--" + INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME, "0"}},
                {new String[]{"--" + INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME, "-1"}},
                {new String[]{"--" + MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_LONG_NAME, "-1"}},
                {new String[]{"--" + NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_PER_FIT_LONG_NAME, "-1"}}
        };
    }

    @Test(dataProvider = "dataACNVArgumentValidation", expectedExceptions = IllegalArgumentException.class)
    public void testACNVArgumentValidation(final String[] testArguments) {
        final String[] commonArguments = {
                "--" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME, COVERAGES_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME, SEGMENT_FILE.getAbsolutePath(),
                "--" + OUTPUT_PREFIX_LONG_NAME, DUMMY_OUTPUT_PREFIX
        };
        final String[] arguments = ArrayUtils.addAll(commonArguments, testArguments);
        testACNV(arguments, DUMMY_OUTPUT_PREFIX);
    }

    @Test
    public void testACNVWithoutAllelicPoN() {
        final File tempDir = createTempDir("acnv-integration-without-pon-" + SAMPLE_NAME);
        final String tempDirPath = tempDir.getAbsolutePath();
        final String outputPrefix = tempDirPath + "/" + SAMPLE_NAME;

        final String[] arguments = {
                "--" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME, COVERAGES_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME, SEGMENT_FILE.getAbsolutePath(),
                "--" + OUTPUT_PREFIX_LONG_NAME, outputPrefix,
                "--" + NUM_SAMPLES_COPY_RATIO_LONG_NAME, "25",
                "--" + NUM_BURN_IN_COPY_RATIO_LONG_NAME, "10",
                "--" + NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME, "25",
                "--" + NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME, "10",
                "--verbosity", "INFO"
        };
        testACNV(arguments, outputPrefix);
    }

    @Test
    public void testACNVWithAllelicPoN() {
        final File tempDir = createTempDir("acnv-integration-with-pon-" + SAMPLE_NAME);
        final String tempDirPath = tempDir.getAbsolutePath();
        final String outputPrefix = tempDirPath + "/" + SAMPLE_NAME;

        final String[] arguments = {
                "--" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME, COVERAGES_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME, SEGMENT_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.ALLELIC_PON_FILE_LONG_NAME, ALLELIC_PON_FILE.getAbsolutePath(),
                "--" + OUTPUT_PREFIX_LONG_NAME, outputPrefix,
                "--" + NUM_SAMPLES_COPY_RATIO_LONG_NAME, "25",
                "--" + NUM_BURN_IN_COPY_RATIO_LONG_NAME, "10",
                "--" + NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME, "25",
                "--" + NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME, "10",
                "--verbosity", "INFO"
        };
        testACNV(arguments, outputPrefix);
    }

    @Test
    public void testACNVWithoutSimilarSegmentMerging() {
        final File tempDir = createTempDir("acnv-integration-without-sim-seg-" + SAMPLE_NAME);
        final String tempDirPath = tempDir.getAbsolutePath();
        final String outputPrefix = tempDirPath + "/" + SAMPLE_NAME;

        final String[] arguments = {
                "--" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME, COVERAGES_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME, SEGMENT_FILE.getAbsolutePath(),
                "--" + OUTPUT_PREFIX_LONG_NAME, outputPrefix,
                "--" + USE_ALL_COPY_RATIO_SEGMENTS_LONG_NAME, "true",
                "--" + NUM_SAMPLES_COPY_RATIO_LONG_NAME, "25",
                "--" + NUM_BURN_IN_COPY_RATIO_LONG_NAME, "10",
                "--" + NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME, "25",
                "--" + NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME, "10",
                "--" + MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_LONG_NAME, "0",
                "--verbosity", "INFO"
        };
        testACNV(arguments, outputPrefix);

        //check that no similar segments were merged
        final File initialSimilarSegmentsFile = new File(outputPrefix + "-" + INITIAL_FIT_FILE_TAG + SEGMENT_FILE_SUFFIX);
        final File finalSimilarSegmentsFile = new File(outputPrefix + "-" + FINAL_FIT_FILE_TAG + SEGMENT_FILE_SUFFIX);
        final List<ACNVModeledSegment> initialACNVModeledSegments = SegmentUtils.readACNVModeledSegmentFile(initialSimilarSegmentsFile);
        final List<ACNVModeledSegment> finalACNVModeledSegments = SegmentUtils.readACNVModeledSegmentFile(finalSimilarSegmentsFile);
        Assert.assertEquals(initialACNVModeledSegments.size(), finalACNVModeledSegments.size());
    }

    @Test
    public void testACNVWithSomeSimilarSegmentMergingRefitting() {
        final File tempDir = createTempDir("acnv-integration-with-some-sim-seg-refitting" + SAMPLE_NAME);
        final String tempDirPath = tempDir.getAbsolutePath();
        final String outputPrefix = tempDirPath + "/" + SAMPLE_NAME;

        final String[] arguments = {
                "--" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME, TUMOR_ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME, COVERAGES_FILE.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME, SEGMENT_FILE.getAbsolutePath(),
                "--" + OUTPUT_PREFIX_LONG_NAME, outputPrefix,
                "--" + USE_ALL_COPY_RATIO_SEGMENTS_LONG_NAME, "true",
                "--" + NUM_SAMPLES_COPY_RATIO_LONG_NAME, "25",
                "--" + NUM_BURN_IN_COPY_RATIO_LONG_NAME, "10",
                "--" + NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME, "25",
                "--" + NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME, "10",
                "--" + NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_PER_FIT_LONG_NAME, "2",
                "--verbosity", "INFO"
        };
        testACNV(arguments, outputPrefix);

        //check that some similar segments were merged
        final File initialSimilarSegmentsFile = new File(outputPrefix + "-" + INITIAL_FIT_FILE_TAG + SEGMENT_FILE_SUFFIX);
        final File finalSimilarSegmentsFile = new File(outputPrefix + "-" + FINAL_FIT_FILE_TAG + SEGMENT_FILE_SUFFIX);
        final List<ACNVModeledSegment> initialACNVModeledSegments = SegmentUtils.readACNVModeledSegmentFile(initialSimilarSegmentsFile);
        final List<ACNVModeledSegment> finalACNVModeledSegments = SegmentUtils.readACNVModeledSegmentFile(finalSimilarSegmentsFile);
        Assert.assertTrue(initialACNVModeledSegments.size() > finalACNVModeledSegments.size());
    }

    private void testACNV(final String[] arguments, final String outputPrefix) {
        runCommandLine(arguments);

        //only check that files are created, do not check for correctness of results
        final File finalSNPSegmentsFile = new File(outputPrefix + "-" + SNP_MAF_SEG_FILE_TAG + SEGMENT_FILE_SUFFIX);
        final File unionedSegmentsFile = new File(outputPrefix + "-" + UNION_SEG_FILE_TAG + SEGMENT_FILE_SUFFIX);
        final File initialSimilarSegmentsFile = new File(outputPrefix + "-" + INITIAL_FIT_FILE_TAG + SEGMENT_FILE_SUFFIX);
        final File finalSimilarSegmentsFile = new File(outputPrefix + "-" + FINAL_FIT_FILE_TAG + SEGMENT_FILE_SUFFIX);
        final File finalCopyRatioParametersFile = new File(outputPrefix + "-" + FINAL_FIT_FILE_TAG + CR_PARAMETER_FILE_SUFFIX);
        final File finalAlleleFractionParametersFile = new File(outputPrefix + "-" + FINAL_FIT_FILE_TAG + AF_PARAMETER_FILE_SUFFIX);

        for (final File outputFile : new File[] {finalSNPSegmentsFile, unionedSegmentsFile,
                initialSimilarSegmentsFile, finalSimilarSegmentsFile, finalCopyRatioParametersFile, finalAlleleFractionParametersFile}) {
            //check that all files are files with a size greater than 0.
            Assert.assertTrue(outputFile.isFile(), outputFile.getAbsolutePath() + " is not a file.");
            Assert.assertTrue(outputFile.length() > 0);
            try {
                //check that all files have:
                //  - at least two lines and all either start with "#" or contain at least one "\t"
                //  - at least two lines with tab (column names + 1 segment or parameter)
                final List<String> outputLines = FileUtils.readLines(outputFile);
                Assert.assertTrue(outputLines.size() >= 2);
                Assert.assertEquals(outputLines.stream().filter(l -> l.contains("\t") || l.startsWith("#")).count(), outputLines.size());
                Assert.assertTrue(outputLines.stream().filter(l -> l.split("\t").length > 2 && !l.startsWith("#")).count() > 2, "File: " + outputFile + " does not seem to have at least one segment and a header.");
                //check that sample names, if present in file, do not contain "/"
                try {
                    final List<String> sampleNamesInOutputFile = SampleNameFinder.determineSampleNamesFromSegmentFile(outputFile);
                    sampleNamesInOutputFile.stream().forEach(n -> Assert.assertFalse(n.contains("/")));
                } catch (final UserException.BadInput e) {
                    //this exception is expected if file does not contain sample name (e.g., ACS seg file), so we do nothing
                }
            } catch (final IOException ioe) {
                Assert.fail("Could not read file: " + outputFile, ioe);
            }
        }
    }
}