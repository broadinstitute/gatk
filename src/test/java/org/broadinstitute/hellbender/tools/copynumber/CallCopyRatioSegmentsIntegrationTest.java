package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Integration test for {@link CallCopyRatioSegments}.
 */
public final class CallCopyRatioSegmentsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber");
    private static final File TEST_SEGMENTS = new File(TEST_SUB_DIR, "call-copy-ratio-segments-segments.seg");

    @Test
    public void testCallSegments() {
        final File outputFile = createTempFile("test.called",".seg");
        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder()
                .addInput(TEST_SEGMENTS)
                .addOutput(outputFile);
        runCommandLine(argsBuilder);
        final CopyRatioSegmentCollection copyRatioSegments = new CopyRatioSegmentCollection(TEST_SEGMENTS);

        final CalledCopyRatioSegmentCollection calledCopyRatioSegments = new CalledCopyRatioSegmentCollection(outputFile);
        Assert.assertEquals(calledCopyRatioSegments.getMetadata(), copyRatioSegments.getMetadata());
        Assert.assertEquals(calledCopyRatioSegments.getIntervals(), copyRatioSegments.getIntervals());
        Assert.assertEquals(calledCopyRatioSegments.getRecords().stream().map(s -> s.getCall().getOutputString()).toArray(), new String[] {"+", "-", "0", "0"});
    }
}
