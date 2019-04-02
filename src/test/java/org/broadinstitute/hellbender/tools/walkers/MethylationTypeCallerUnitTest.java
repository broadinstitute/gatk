package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.List;

public class MethylationTypeCallerUnitTest extends CommandLineProgramTest {

    // If true, update the expected outputs in tests that assert an exact match vs. prior output,
    // instead of actually running the tests. Can be used with "./gradlew test -Dtest.single=MethylationTypeCallerUnitTest"
    // to update all of the exact-match tests at once. After you do this, you should look at the
    // diffs in the new expected outputs in git to confirm that they are consistent with expectations.
    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    public static final String TEST_FILES_DIR = toolsTestDir + "MethylationTypeCaller/";
    public static final String TEST_FILES_INPUT_REFERENCE_DIR = largeFileTestDir + "GRCm38_primary_assembly_genome/";

    /*
     * Test that in VCF mode we're consistent with past  results
     */
    @Test(dataProvider="getMethylationTypeCallerTestInput")
    public void testVCFIsConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File outputVCF = createTempFile("testVCFIsConsistentWithPastResults", ".vcf");
        final File expectedVCF = new File(TEST_FILES_DIR, "chr14.unique_reads.methylC_seq.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedVCF.getAbsolutePath() : outputVCF.getAbsolutePath();

        final String[] inputArgs = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-O", outputPath
        };

        runCommandLine(inputArgs);

        // Test for an exact match against past results
        if (!UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            final List<VariantContext> outputVariants = VariantContextTestUtils.readEntireVCFIntoMemory(outputVCF.toString()).getRight();
            final List<VariantContext> expectedVariants = VariantContextTestUtils.readEntireVCFIntoMemory(expectedVCF.toString()).getRight();
            Assert.assertEquals(expectedVariants.toString(), outputVariants.toString());
        }

    }

    @DataProvider
    public Object[][] getMethylationTypeCallerTestInput() {
        final String inputFileName = TEST_FILES_DIR + "chr14.unique_reads.methylC_seq.bam";
        final String referenceFileName = TEST_FILES_INPUT_REFERENCE_DIR + "chr14.GRCm38.primary_assembly.genome.fa.gz";
        return new Object[][] {
                {inputFileName, referenceFileName}
        };
    }
}

