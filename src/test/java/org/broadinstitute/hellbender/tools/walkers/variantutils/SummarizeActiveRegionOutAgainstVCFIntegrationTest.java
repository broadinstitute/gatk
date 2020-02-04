package org.broadinstitute.hellbender.tools.walkers.variantutils;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;

public class SummarizeActiveRegionOutAgainstVCFIntegrationTest extends CommandLineProgramTest {
    public static final String TEST_FILES_DIR = toolsTestDir + "haplotypecaller/";

    @Test
    // NOTE these results should be updated in the future if there are massive changes to the output of the haplotype caller outputs
    public void testAnnotationOfAssemblyRegion() throws Exception {
        final File expected = new File(TEST_FILES_DIR, "testAnnotationOfAssemblyRegionOut.expected.table");
        final File output = createTempFile("testAnnotationOfAssemblyRegion", ".table");
        final File inputAssemblyRegion = new File(TEST_FILES_DIR, "testAssemblyRegionOutput.table");
        final File inputVCF = new File(TEST_FILES_DIR, "expected.testVCFMode.gatk4.vcf");

        final String[] args = {
                "--active-region-summary", inputAssemblyRegion.getAbsolutePath(),
                "--variants-to-count", inputVCF.getAbsolutePath(),
                "--output", output.getAbsolutePath()
        };

        runCommandLine(args);

        // Test for an exact match against past results
        IntegrationTestSpec.assertEqualTextFiles(output, expected);

    }

}