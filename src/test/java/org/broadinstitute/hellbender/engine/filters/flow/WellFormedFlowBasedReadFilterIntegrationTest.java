package org.broadinstitute.hellbender.engine.filters.flow;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.walkers.variantrecalling.FlowTestConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class WellFormedFlowBasedReadFilterIntegrationTest extends CommandLineProgramTest {

    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    private static String testDir = publicTestDir + FlowTestConstants.READ_FILTER_DATA_DIR;

    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @Test
    public void testBasic() throws IOException {

        final File outputDir = createTempDir("testReadFilterTest");
        final File expectedFile = new File(testDir + "/read_filter_output.sam");
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/read_filter_output.sam");
        final File input = new File(largeFileTestDir, "input_jukebox_for_test.bam");

        final String[] args = new String[]{
                "-I", input.getAbsolutePath(),
                "-O", outputFile.getAbsolutePath(),
                "--intervals", "chr9:81149486-81177047",
                "--read-filter", WellformedFlowBasedReadFilter.class.getSimpleName()
        };

        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the otuput file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        if (!UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile, "@");
        }
    }

    @Override
    public String getTestedToolName() {
        return "PrintReads";
    }
}
