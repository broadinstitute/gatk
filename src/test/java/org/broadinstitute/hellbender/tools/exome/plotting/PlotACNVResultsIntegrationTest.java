package org.broadinstitute.hellbender.tools.exome.plotting;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class PlotACNVResultsIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";

    //test files
    private static final File SNP_COUNTS_FILE = new File(TEST_SUB_DIR, "snps-for-plotting.tsv");
    private static final File TANGENT_NORMALIZED_COUNTS_FILE = new File(TEST_SUB_DIR, "coverages-for-plotting.tsv");
    private static final File SEGMENTS_FILE = new File(TEST_SUB_DIR, "acnv-segments-for-plotting.seg");
    private static final File SEQUENCE_DICTIONARY_FILE = new File(TEST_SUB_DIR, "sequence-dictionary-for-plotting.dict");

    //test files for invalid configurations
    private static final File SEQUENCE_DICTIONARY_NO_CONTIGS_ABOVE_MINIMUM_LENGTH_FILE = new File(TEST_SUB_DIR, "sequence-dictionary-for-plotting-no-contigs-above-minimum-length.dict");
    private static final File COUNTS_BAD_SAMPLE_NAME_FILE = new File(TEST_SUB_DIR, "coverages-for-plotting-bad-sample-name.tsv");
    private static final File SEGMENTS_BAD_SAMPLE_NAME_FILE = new File(TEST_SUB_DIR, "acnv-segments-for-plotting-bad-sample-name.seg");
    private static final File SNP_COUNTS_DATA_OUT_OF_BOUNDS_FILE = new File(TEST_SUB_DIR, "snps-for-plotting-data-out-of-bounds.tsv");
    private static final File COUNTS_DATA_OUT_OF_BOUNDS_FILE = new File(TEST_SUB_DIR, "coverages-for-plotting-data-out-of-bounds.tsv");
    private static final File SEGMENTS_DATA_OUT_OF_BOUNDS_FILE = new File(TEST_SUB_DIR, "acnv-segments-for-plotting-data-out-of-bounds.seg");
    
    private static final String OUTPUT_PREFIX = "test";
    private static final int THRESHOLD_PLOT_FILE_SIZE_IN_BYTES = 50000;  //test that data points are plotted (not just background/axes)

    //checks that output files with reasonable file sizes are generated, but correctness of output is not checked
    @Test
    public void testACNVPlotting() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_SHORT_NAME, SNP_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, TANGENT_NORMALIZED_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEGMENTS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + PlotACNVResults.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX,
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + "_ACNV.png").exists());
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + "_ACNV.png").length() > THRESHOLD_PLOT_FILE_SIZE_IN_BYTES);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMinimumContigLength() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_SHORT_NAME, SNP_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, TANGENT_NORMALIZED_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEGMENTS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_NO_CONTIGS_ABOVE_MINIMUM_LENGTH_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + PlotACNVResults.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.class)
    public void testOutputDirExists() throws IOException {
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_SHORT_NAME, SNP_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, TANGENT_NORMALIZED_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEGMENTS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "Non-existent-path",
                "-" + PlotACNVResults.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.class)
    public void testMissingSNPCountsFile() throws IOException {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_SHORT_NAME, "Non-existent-file",
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, TANGENT_NORMALIZED_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEGMENTS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + PlotACNVResults.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.class)
    public void testMissingTangentFile() throws IOException {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_SHORT_NAME, SNP_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,  "Non-existent-file",
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEGMENTS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + PlotACNVResults.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.class)
    public void testMissingSegmentsFile() throws IOException {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_SHORT_NAME, SNP_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, TANGENT_NORMALIZED_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,  "Non-existent-file",
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + PlotACNVResults.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.class)
    public void testMissingSequenceDictionaryFile() throws IOException {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_SHORT_NAME, SNP_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, TANGENT_NORMALIZED_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEGMENTS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME,  "Non-existent-file",
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + PlotACNVResults.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testTangentSampleNameMismatch() throws IOException {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_SHORT_NAME, SNP_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, COUNTS_BAD_SAMPLE_NAME_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEGMENTS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + PlotSegmentedCopyRatio.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testSegmentsSampleNameMismatch() throws IOException {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_SHORT_NAME, SNP_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, TANGENT_NORMALIZED_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEGMENTS_BAD_SAMPLE_NAME_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + PlotSegmentedCopyRatio.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testSNPCountsDataOutOfBounds() throws IOException {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_SHORT_NAME, SNP_COUNTS_DATA_OUT_OF_BOUNDS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, TANGENT_NORMALIZED_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEGMENTS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + PlotSegmentedCopyRatio.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testTangentDataOutOfBounds() throws IOException {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_SHORT_NAME, SNP_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, COUNTS_DATA_OUT_OF_BOUNDS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEGMENTS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + PlotSegmentedCopyRatio.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testSegmentsDataOutOfBounds() throws IOException {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.ALLELIC_COUNTS_FILE_SHORT_NAME, SNP_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, TANGENT_NORMALIZED_COUNTS_FILE.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, SEGMENTS_DATA_OUT_OF_BOUNDS_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "-" + PlotSegmentedCopyRatio.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }
}