package org.broadinstitute.hellbender.tools.copynumber.plotting;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class PlotModeledSegmentsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber/plotting");

    //test files
    private static final File DENOISED_COPY_RATIOS_FILE = new File(TEST_SUB_DIR, "plotting-copy-ratios.tsv");
    private static final File ALLELIC_COUNTS_FILE = new File(TEST_SUB_DIR, "plotting-allelic-counts.tsv");
    private static final File MODELED_SEGMENTS_FILE = new File(TEST_SUB_DIR, "plotting-modeled-segments.seg");
    private static final File SEQUENCE_DICTIONARY_FILE = new File(TEST_SUB_DIR, "plotting-sequence-dictionary.dict");

    //test files for invalid configurations
    private static final File SEQUENCE_DICTIONARY_WITH_NO_CONTIGS_ABOVE_MINIMUM_LENGTH_FILE = new File(TEST_SUB_DIR, "plotting-sequence-dictionary-with-no-contigs-above-minimum-length.dict");
    private static final File DENOISED_COPY_RATIOS_WITH_SAMPLE_NAME_MISMATCH_FILE = new File(TEST_SUB_DIR, "plotting-copy-ratios-with-sample-name-mismatch.tsv");
    private static final File DENOISED_COPY_RATIOS_OUT_OF_DICTIONARY_BOUNDS_FILE = new File(TEST_SUB_DIR, "plotting-copy-ratios-out-of-dictionary-bounds.tsv");
    private static final File ALLELIC_COUNTS_WITH_SAMPLE_NAME_MISMATCH_FILE = new File(TEST_SUB_DIR, "plotting-allelic-counts-with-sample-name-mismatch.tsv");
    private static final File ALLELIC_COUNTS_OUT_OF_DICTIONARY_BOUNDS_FILE = new File(TEST_SUB_DIR, "plotting-allelic-counts-out-of-dictionary-bounds.tsv");
    private static final File MODELED_SEGMENTS_WITH_SAMPLE_NAME_MISMATCH_FILE = new File(TEST_SUB_DIR, "plotting-modeled-segments-with-sample-name-mismatch.seg");
    private static final File MODELED_SEGMENTS_OUT_OF_DICTIONARY_BOUNDS_FILE = new File(TEST_SUB_DIR, "plotting-modeled-segments-out-of-dictionary-bounds.seg");
    private static final File MODELED_SEGMENTS_WITH_WRONG_NUM_POINTS_COPY_RATIO_FILE = new File(TEST_SUB_DIR, "plotting-modeled-segments-with-wrong-num-points-copy-ratio.seg");
    private static final File MODELED_SEGMENTS_WITH_WRONG_NUM_POINTS_ALLELE_FRACTION_FILE = new File(TEST_SUB_DIR, "plotting-modeled-segments-with-wrong-num-points-allele-fraction.seg");
    
    private static final String OUTPUT_PREFIX = "test";
    private static final int THRESHOLD_PLOT_FILE_SIZE_IN_BYTES = 50000;  //test that data points are plotted (not just background/axes)

    //checks that output files with reasonable file sizes are generated, but correctness of output is not checked
    @Test
    public void testPlotting() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "--" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME, MODELED_SEGMENTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + ".modeled.png").exists());
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + ".modeled.png").length() > THRESHOLD_PLOT_FILE_SIZE_IN_BYTES);
    }

    @Test
    public void testPlottingDenoisedCopyRatiosOnly() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "--" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME, MODELED_SEGMENTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + ".modeled.png").exists());
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + ".modeled.png").length() > THRESHOLD_PLOT_FILE_SIZE_IN_BYTES / 2);    //copy-ratio-only plot is half the size
    }

    @Test
    public void testPlottingAllelicCountsOnly() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "--" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME, MODELED_SEGMENTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + ".modeled.png").exists());
        Assert.assertTrue(new File(outputDir, OUTPUT_PREFIX + ".modeled.png").length() > THRESHOLD_PLOT_FILE_SIZE_IN_BYTES / 2);    //allele-fraction-only plot is half the size
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMinimumContigLength() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "--" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME, MODELED_SEGMENTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, SEQUENCE_DICTIONARY_WITH_NO_CONTIGS_ABOVE_MINIMUM_LENGTH_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.class)
    public void testOutputDirExists() {
        final String[] arguments = {
                "--" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME, MODELED_SEGMENTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "Non-existent-path",
                "--" + CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.class)
    public void testNonExistentDenoisedCopyRatiosFile() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "--" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, "Non-existent-file",
                "--" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME, MODELED_SEGMENTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.class)
    public void testNonExistentAllelicCountsFile() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "--" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, "Non-existent-file",
                "--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME, MODELED_SEGMENTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.class)
    public void testNonExistentModeledSegmentsFile() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "--" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME, "Non-existent-file",
                "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = UserException.class)
    public void testNonExistentSequenceDictionaryFile() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "--" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME, MODELED_SEGMENTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, "Non-existent-file",
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDenoisedCopyRatiosSampleNameMismatch() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "--" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, DENOISED_COPY_RATIOS_WITH_SAMPLE_NAME_MISMATCH_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME, MODELED_SEGMENTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAllelicCountsSampleNameMismatch() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "--" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, ALLELIC_COUNTS_WITH_SAMPLE_NAME_MISMATCH_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME, MODELED_SEGMENTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testModeledSegmentsSampleNameMismatch() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "--" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME, MODELED_SEGMENTS_WITH_SAMPLE_NAME_MISMATCH_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDenoisedCopyRatiosOutOfDictionaryBounds() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "--" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, DENOISED_COPY_RATIOS_OUT_OF_DICTIONARY_BOUNDS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME, MODELED_SEGMENTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAllelicCountsOutOfDictionaryBounds() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "--" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, ALLELIC_COUNTS_OUT_OF_DICTIONARY_BOUNDS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME, MODELED_SEGMENTS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testModeledSegmentsOutOfDictionaryBounds() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "--" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME, MODELED_SEGMENTS_OUT_OF_DICTIONARY_BOUNDS_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testModeledSegmentsWithWrongNumPointsCopyRatio() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "--" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME, MODELED_SEGMENTS_WITH_WRONG_NUM_POINTS_COPY_RATIO_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testModeledSegmentsWithWrongNumPointsAlleleFraction() {
        final File outputDir = createTempDir("testDir");
        final String[] arguments = {
                "--" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME, DENOISED_COPY_RATIOS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME, ALLELIC_COUNTS_FILE.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME, MODELED_SEGMENTS_WITH_WRONG_NUM_POINTS_ALLELE_FRACTION_FILE.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, SEQUENCE_DICTIONARY_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputDir.getAbsolutePath(),
                "--" + CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME, OUTPUT_PREFIX
        };
        runCommandLine(arguments);
    }
}