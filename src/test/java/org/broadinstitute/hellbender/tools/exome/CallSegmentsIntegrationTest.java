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
    private static final File TEST_TARGETS = new File(TEST_DIR, "targets.tsv");
    private static final File TEST_SEGMENTS = new File(TEST_DIR, "segments.tsv");
    private static final File TEST_SEGMENTS_LEGACY = new File(TEST_DIR, "segments_legacy.tsv");

    @Test
    public void testCallSegments() {
        final File outputFile = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, TEST_SEGMENTS.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, TEST_TARGETS.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath()
        };
        runCommandLine(arguments);

        final List<ModeledSegment> calls = SegmentUtils.readModeledSegmentsFromSegmentFile(outputFile);
        Assert.assertEquals(calls.stream().map(ModeledSegment::getCall).toArray(), new String[] {"+", "-", "0", "0"});
    }

    @Test
    public void testCallSegmentsLegacy() {
        final File outputFile = createTempFile("test",".txt");

        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, TEST_SEGMENTS_LEGACY.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME, TEST_TARGETS.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + ExomeStandardArgumentDefinitions.LEGACY_SEG_FILE_SHORT_NAME
        };
        runCommandLine(arguments);

        final List<ModeledSegment> calls = SegmentUtils.readModeledSegmentsFromSegmentFile(outputFile);
        Assert.assertEquals(calls.stream().map(ModeledSegment::getCall).toArray(), new String[] {"+", "-", "0", "0"});
    }
}
