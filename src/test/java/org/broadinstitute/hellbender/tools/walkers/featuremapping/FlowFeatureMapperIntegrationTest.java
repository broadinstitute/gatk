package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.tools.walkers.variantrecalling.FlowTestConstants;

import java.io.File;
import java.io.IOException;

public class FlowFeatureMapperIntegrationTest extends CommandLineProgramTest {

    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    private static String testDir = publicTestDir + FlowTestConstants.FEATURE_MAPPING_DATA_DIR;

    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @Test
    public void testBasic() throws IOException {

        final File outputDir = createTempDir("testFlowFeatureMapperTest");
        final File expectedFile = new File(testDir + "/snv_feature_mapper_output.vcf");
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/snv_feature_mapper_output.vcf");

        final String[] args = new String[] {
                "-R", largeFileTestDir + "/Homo_sapiens_assembly38.fasta.gz",
                "-O", outputFile.getAbsolutePath(),
                "-I", testDir + "/snv_feature_mapper_input.bam",
                "--copy-attr", "tr",
                "--limit-score", "100",
                "--min-score", "0",
                "--snv-identical-bases", "10",
                "--debug-negatives", "false",
                "--debug-read-name", "150451-BC94-0645901755"
        };

        // run the tool
        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the otuput file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile, "#");
        }
    }

}
