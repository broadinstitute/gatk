package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

/**
 * Integration test for {@link CallSegments}.
 */
public final class CallSegmentsIntegrationTest extends CommandLineProgramTest{
    private static final File TEST_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome/caller");
    private static final File TEST_TARGETS = new File(TEST_DIR,"targets.tsv");
    private static final File TEST_SEGMENTS = new File(TEST_DIR,"segments.tsv");
    private static final File TEST_SEGMENTS_LEGACY = new File(TEST_DIR,"segments_legacy.tsv");
    private static final String SAMPLE_NAME = "sample";

    @Test
    public void testCallSegmentsExperimental() {
        final File outputFile = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, TEST_SEGMENTS.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, TEST_TARGETS.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.SAMPLE_LONG_NAME, SAMPLE_NAME,
                "-" + CallSegments.Z_THRESHOLD_SHORT_NAME, Double.toString(CnvCaller.DEFAULT_Z_SCORE_THRESHOLD)
        };
        runCommandLine(arguments);

        final List<ModeledSegment> calls = SegmentUtils.readModeledSegmentsFromSegmentFile(outputFile);

        Assert.assertEquals(calls.get(0).getCall(), "+");
        Assert.assertEquals(calls.get(1).getCall(), "-");
        Assert.assertEquals(calls.get(2).getCall(), "0");
        Assert.assertEquals(calls.get(3).getCall(), "0");
    }

    @Test
    public void testCallSegmentsLegacyExperimental() {
        final File outputFile = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, TEST_SEGMENTS_LEGACY.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, TEST_TARGETS.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.SAMPLE_LONG_NAME, SAMPLE_NAME,
                "-" + CallSegments.Z_THRESHOLD_SHORT_NAME, Double.toString(CnvCaller.DEFAULT_Z_SCORE_THRESHOLD),
                "-" + ExomeStandardArgumentDefinitions.LEGACY_SEG_FILE_SHORT_NAME, "-" + CallSegments.EXPERIMENTAL_CALLER_SHORT_NAME
        };
        runCommandLine(arguments);

        final List<ModeledSegment> calls = SegmentUtils.readModeledSegmentsFromSegmentFile(outputFile);

        Assert.assertEquals(calls.get(0).getCall(), "+");
        Assert.assertEquals(calls.get(1).getCall(), "-");
        Assert.assertEquals(calls.get(2).getCall(), "0");
        Assert.assertEquals(calls.get(3).getCall(), "0");
    }

    @Test
    public void testCallSegments() {
        final File outputFile = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, TEST_SEGMENTS.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, TEST_TARGETS.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.SAMPLE_LONG_NAME, SAMPLE_NAME,
                "-" + CallSegments.Z_THRESHOLD_SHORT_NAME, Double.toString(CnvCaller.DEFAULT_Z_SCORE_THRESHOLD),
                "-" + CallSegments.EXPERIMENTAL_CALLER_SHORT_NAME
        };
        runCommandLine(arguments);

        final List<ModeledSegment> calls = SegmentUtils.readModeledSegmentsFromSegmentFile(outputFile);

        Assert.assertEquals(calls.get(0).getCall(), "+");
        Assert.assertEquals(calls.get(1).getCall(), "-");
        Assert.assertEquals(calls.get(2).getCall(), "0");
        Assert.assertEquals(calls.get(3).getCall(), "0");
    }

    @Test
    public void testCallSegmentsLegacy() {
        final File outputFile = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, TEST_SEGMENTS_LEGACY.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME, TEST_TARGETS.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "--" + ExomeStandardArgumentDefinitions.SAMPLE_LONG_NAME, SAMPLE_NAME,
                "-" + CallSegments.Z_THRESHOLD_SHORT_NAME, Double.toString(CnvCaller.DEFAULT_Z_SCORE_THRESHOLD),
                "-" + ExomeStandardArgumentDefinitions.LEGACY_SEG_FILE_SHORT_NAME
        };
        runCommandLine(arguments);

        final List<ModeledSegment> calls = SegmentUtils.readModeledSegmentsFromSegmentFile(outputFile);

        Assert.assertEquals(calls.get(0).getCall(), "+");
        Assert.assertEquals(calls.get(1).getCall(), "-");
        Assert.assertEquals(calls.get(2).getCall(), "0");
        Assert.assertEquals(calls.get(3).getCall(), "0");
    }
}
