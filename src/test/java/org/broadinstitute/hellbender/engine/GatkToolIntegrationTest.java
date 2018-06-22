package org.broadinstitute.hellbender.engine;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class GatkToolIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";

    @Test
    public void testSitesOnlyMode() {
        File out = createTempFile("GTStrippedOutput", "vcf");
        String[] args = new String[] {
                "-V",  TEST_DIRECTORY + "vcf_with_genotypes.vcf",
                "--" + StandardArgumentDefinitions.SITES_ONLY_LONG_NAME,
                "-O",
                out.getAbsolutePath()};
        runCommandLine(Arrays.asList(args), SelectVariants.class.getSimpleName());

        // Assert that the genotype field has been stripped from the file
        Pair<VCFHeader, List<VariantContext>> results = VariantContextTestUtils.readEntireVCFIntoMemory(out.getAbsolutePath());

        Assert.assertFalse(results.getLeft().hasGenotypingData());
        for (VariantContext v: results.getRight()) {
            Assert.assertFalse(v.hasGenotypes());
        }
    }
}
