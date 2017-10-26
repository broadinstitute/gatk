package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.coverage.caller.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberStandardArgument;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Integration test for {@link CallCopyRatioSegments}.
 */
public final class CallCopyRatioSegmentsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DIR = new File(toolsTestDir, "copynumber/coverage/caller");
    private static final File TEST_DENOISED_COPY_RATIOS = new File(TEST_DIR, "call-copy-ratio-segments-denoised-copy-ratios.tsv");
    private static final File TEST_SEGMENTS = new File(TEST_DIR, "call-copy-ratio-segments-segments.seg");

    @Test
    public void testCallSegments() {
        final File outputFile = createTempFile("test.called",".seg");

        final String[] arguments = {
                "-" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, TEST_DENOISED_COPY_RATIOS.getAbsolutePath(),
                "-" + CopyNumberStandardArgument.SEGMENTS_FILE_SHORT_NAME, TEST_SEGMENTS.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath()
        };
        runCommandLine(arguments);

        final CalledCopyRatioSegmentCollection calledCopyRatioSegments = new CalledCopyRatioSegmentCollection(outputFile);
        Assert.assertEquals(calledCopyRatioSegments.getRecords().stream().map(s -> s.getCall().getOutputString()).toArray(), new String[] {"+", "-", "0", "0"});
    }
}
