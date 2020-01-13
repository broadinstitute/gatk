package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public class MTLowHeteroplasmyFilterTest extends CommandLineProgramTest {
    private static final File MITO_REF = new File(toolsTestDir, "mutect/mito/Homo_sapiens_assembly38.mt_only.fasta");
    private static final File NA12878_MITO_FILTERED_VCF = new File(toolsTestDir, "mutect/mito/filtered.vcf");

    @Test
    public void testLowHetVariantWalker() throws IOException {
        final IntegrationTestSpec testSpec = new IntegrationTestSpec(
                        " -R " + MITO_REF.getAbsolutePath()  +
                        " -V " + NA12878_MITO_FILTERED_VCF +
                        " -O %s",
                Arrays.asList(toolsTestDir + "mutect/mito/expected_LowHetVariantWalkerIntegrationTest_output.txt")
        );
        testSpec.executeTest("testLowHetVariantWalker", this);

        final IntegrationTestSpec testLowHetNoneSpec = new IntegrationTestSpec(
                " -R " + MITO_REF.getAbsolutePath()  +
                        " -V " + NA12878_MITO_FILTERED_VCF +
                        " -O %s" +
                        " --min-low-het-sites 5",
                Arrays.asList(toolsTestDir + "mutect/mito/expected_LowHetNone_output.txt")
        );
        testLowHetNoneSpec.executeTest("testLowHetVariantWalker", this);
    }

}
