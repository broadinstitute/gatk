package org.broadinstitute.hellbender.tools.walkers.annotator;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GatkToolIntegrationTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.mutect.SomaticLikelihoodsEngine;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

/**
 * Integration tests for the {@link AllelePseudoDepth} annotator.
 */
public class DirichletAlleleDepthAndFractionIntegrationTest extends GatkToolIntegrationTest {

    public static final String TEST_FILES_DIR = toolsTestDir + "haplotypecaller/";

    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = true;

    @Override
    public String getTestedToolName() {
        return HaplotypeCaller.class.getSimpleName();
    }

    /*
     * Test that in VCF mode we're consistent with past GATK4 results
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testVCFModeIsConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConsistentWithPastResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testVCFMode.gatk4.withDDandDF.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-A", AllelePseudoDepth.class.getSimpleName(),
                "-O", outputPath,
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    @Test
    public void testNoInformativeMatrix() {
        final Array2DRowRealMatrix rm = new Array2DRowRealMatrix(2, 100);
        final double[] postriors = SomaticLikelihoodsEngine.alleleFractionsPosterior(rm, new double[] { 0.001, 0.001});
        System.err.println(Arrays.toString(postriors));
    }

    @DataProvider(name="HaplotypeCallerTestInputs")
    public Object[][] getHaplotypCallerTestInputs() {
        return new Object[][] {
                {NA12878_20_21_WGS_bam, b37_reference_20_21},
                {NA12878_20_21_WGS_cram, b37_reference_20_21}
        };
    }
}
