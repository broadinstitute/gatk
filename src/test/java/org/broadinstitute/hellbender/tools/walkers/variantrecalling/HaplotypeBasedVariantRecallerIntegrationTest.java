package org.broadinstitute.hellbender.tools.walkers.variantrecalling;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class HaplotypeBasedVariantRecallerIntegrationTest extends CommandLineProgramTest {

    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;
    private static String testDir = publicTestDir + FlowTestConstants.VARIANT_CALLING_DATA_DIR;

    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @Test
    public void testBasic() throws IOException {

        final File outputDir = createTempDir("testHaplotypeBasedVariantRecallerTest");
        final File expectedFile = new File(testDir + "/variantRecallerBasic.expected.csv");
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/output.csv");

        final String[] args = new String[] {
                "--alleles-file-vcf", testDir + "/150292-BC05.vcf.gz",
                "--haplotypes-file-bam", testDir + "/haps_chr5.bam",
                "-I", testDir + "/chr5.bam1.rename.bam",
                "-I", testDir + "/chr5.bam2.rename.bam",
                "--reference", largeFileTestDir + "/Homo_sapiens_assembly38.fasta.gz",
                "--matrix-file-csv", outputFile.getAbsolutePath(),
                "--likelihood-calculation-engine", "FlowBased",
                "--phred-scaled-global-read-mismapping-rate", "-1",
                "-L", "chr5:70036625"
        };

        runCommandLine(args);  // no assert, just make sure we don't throw

        // verify that output file has been created
        // walk the output and expected files, compare non-comment lines
        Assert.assertTrue(outputFile.exists());
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile, "#");
        }
    }
}
