package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;
import java.util.List;

public class RefineComplexVariantsIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_DIR = toolsTestDir + "walkers/sv/RefineComplexVariants/";

    @Test
    public void testRefineComplexVariantsSubsetMatchesExpectedOutput() {
        final File output = createTempFile("refine-complex-variants.", ".vcf.gz");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(new File(TEST_DIR + "input.subset.vcf.gz"))
                .add(RefineComplexVariants.BATCH_SAMPLE_LISTS_LONG_NAME, TEST_DIR + "batch_samples.txt")
                .add(RefineComplexVariants.DISCORDANT_PAIR_FILES_LONG_NAME, TEST_DIR + "evidence.pe.txt.gz")
                .add(RefineComplexVariants.DEPTH_DEL_BEDS_LONG_NAME, TEST_DIR + "depth.DEL.bed")
                .add(RefineComplexVariants.DEPTH_DUP_BEDS_LONG_NAME, TEST_DIR + "depth.DUP.bed")
                .addOutput(output);

        runCommandLine(args, RefineComplexVariants.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> expected =
                VariantContextTestUtils.readEntireVCFIntoMemory(TEST_DIR + "expected.subset.vcf.gz");
        final Pair<VCFHeader, List<VariantContext>> actual =
                VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());

        Assert.assertEquals(actual.getRight().size(), expected.getRight().size(),
                "Output and expected VCFs contain different variant counts");

        for (int i = 0; i < expected.getRight().size(); i++) {
            VariantContextTestUtils.assertVariantContextsAreEqual(
                    actual.getRight().get(i),
                    expected.getRight().get(i),
                    Collections.emptyList(),
                    Collections.emptyList());
        }
    }
}
