package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.cram.build.CramContainerIterator;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

public class SplitCRAMIntegrationTest extends CommandLineProgramTest {

    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    private static String testDir = publicTestDir;

    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @Test
    public void testBasic() throws IOException {

        final File outputDir = createTempDir("testSplitCRAM");
        final int OUTPUT_FILE_COUNT = 2;
        final File[] expectedFiles = new File[OUTPUT_FILE_COUNT];
        final File[] outputFiles = new File[OUTPUT_FILE_COUNT];
        for ( int i = 0 ; i < OUTPUT_FILE_COUNT ; i++ ) {
            expectedFiles[i] = new File(testDir + String.format("/large/expected_SplitCRAM_output_%04d.cram", i));
            outputFiles[i] = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFiles[i] : new File(outputDir + String.format("/expected_SplitCRAM_output_%04d.cram", i));
        }
        final String outputTemplate = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS
                ? (testDir + "/large/expected_SplitCRAM_output_%04d.cram")
                : (outputDir + "/expected_SplitCRAM_output_%04d.cram");


        final String[] args = new String[] {
            "-I", testDir + "/large/cnv_somatic_workflows_test_files/SM-74NEG-v1-chr20-downsampled.deduplicated.cram",
            "-O", outputTemplate,
            "--shard-records", "7000"
        };

        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the output file
        for ( File outputFile : outputFiles ) {
            Assert.assertTrue(outputFile.exists());
            SamAssertionUtils.assertCRAMContents(outputFile.toPath());
        }

        // walk the output and expected files, compare non-comment lines
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {

            for ( int i = 0 ; i < OUTPUT_FILE_COUNT ; i++ ) {

                // verify length
                Assert.assertEquals(outputFiles[i].length(), expectedFiles[i].length());

                // check that files can be consumed as CRAM files
                try (final CramContainerIterator iter = new CramContainerIterator(new BufferedInputStream(new FileInputStream(outputFiles[i])))) {
                    Assert.assertNotNull(iter.getCramHeader());
                    Assert.assertNotNull(iter.getSamFileHeader());
                }
            }
        }
    }
}
