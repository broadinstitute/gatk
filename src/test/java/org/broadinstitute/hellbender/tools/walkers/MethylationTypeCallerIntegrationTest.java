package org.broadinstitute.hellbender.tools.walkers;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class MethylationTypeCallerIntegrationTest extends CommandLineProgramTest {

    public static final String TEST_FILES_INPUT_REFERENCE_GRCM38_DIR = largeFileTestDir + "GRCm38_primary_assembly_genome/";

    @Test(dataProvider="getMethylationTypeCallerTestInput")
    public void testBasicMethylationCoverage(final String inputFileName, final String referenceFileName) throws Exception {

        final File outputVCF = createTempFile("testBasicMethylationCoverageVCF", ".vcf");
        final File expectedVCF = getTestFile("chr14_subset.methylC_seq.expected.vcf");

        final String[] inputArgs = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-O", outputVCF.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(inputArgs);

        // Test for an exact match against past results
        IntegrationTestSpec.assertEqualTextFiles(outputVCF, expectedVCF);
    }

    @DataProvider
    public Object[][] getMethylationTypeCallerTestInput() {
        final String inputFileName = getToolTestDataDir() + "chr14_subset.methylC_seq.bam";
        final String referenceFileName = largeFileTestDir + "GRCm38_primary_assembly_genome/chr14.GRCm38.primary_assembly.genome.fa.gz";
        return new Object[][] {
                {inputFileName, referenceFileName}
        };
    }
}

