package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;
import java.util.Iterator;
import java.util.stream.Collectors;

public class NVScoreVariantsIntegrationTest extends CommandLineProgramTest {
    private static final String inputVCF = largeFileTestDir + "VQSR/recalibrated_chr20_start.vcf";
    private static final String inputBAM = largeFileTestDir + "VQSR/g94982_contig_20_start_bamout.bam";
    private static final String reference = b37_reference_20_21;

    private static final double EPSILON_FOR_1D = 0.01;
    private static final double EPSILON_FOR_2D = 0.5;

    // This test for the 1D model PASSES when run locally in the scripts/nvscorevariants_environment.yml
    // conda environment, but cannot be enabled until that conda environment is incorporated into
    // the main GATK conda environment.
    @Test(groups = {"python"}, enabled = false)
    public void test1DModel() {
        final File tempVcf = createTempFile("test1DModel", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf");

        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getAbsolutePath())
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21);

        runCommandLine(argsBuilder);
        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_1D_KEY, EPSILON_FOR_1D);
    }

    // This test for the 2D model FAILS when run locally in the scripts/nvscorevariants_environment.yml
    // conda environment, despite the much higher epsilon of 0.5. There are scores that mismatch by as
    // much as 5+:
    //
    // java.lang.AssertionError: scores at 20:61098 differed by 5.085, which is greater than the allowed tolerance of 0.5
    //
    // This test also cannot be enabled until the nvscorevariants conda environment is incorporated into
    // the main GATK conda environment.
    @Test(groups = {"python"}, enabled = false)
    public void test2DModel() {
        final File tempVcf = createTempFile("test2DModel", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf");

        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getAbsolutePath())
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, b37_reference_20_21)
                .add(StandardArgumentDefinitions.INPUT_LONG_NAME, inputBAM)
                .add("tensor-type", NVScoreVariants.TensorType.read_tensor.name());

        runCommandLine(argsBuilder);
        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_2D_KEY, EPSILON_FOR_2D);
    }

    private void assertInfoFieldsAreClose(final File actualVcf, final File expectedVcf, final String infoKey, final double epsilon) {
        Iterator<VariantContext> expectedVi = VariantContextTestUtils.streamVcf(expectedVcf).collect(Collectors.toList()).iterator();
        Iterator<VariantContext> actualVi = VariantContextTestUtils.streamVcf(actualVcf).collect(Collectors.toList()).iterator();

        while (expectedVi.hasNext() && actualVi.hasNext()) {
            VariantContext expectedVc = expectedVi.next();
            VariantContext actualVc = actualVi.next();
            Assert.assertEquals(actualVc.getContig(), expectedVc.getContig(), "Variants from actual and expected VCFs do not match in their location");
            Assert.assertEquals(actualVc.getStart(), expectedVc.getStart(), "Variants from actual and expected VCFs do not match in their location");
            double expectedScore = expectedVc.getAttributeAsDouble(infoKey, 0.0); // Different defaults trigger failures on missing scores
            double actualScore = actualVc.getAttributeAsDouble(infoKey, epsilon+1.0);
            double diff = Math.abs(expectedScore-actualScore);
            Assert.assertTrue(diff < epsilon, "scores at " + expectedVc.getContig() + ":" + expectedVc.getStart() +
                    " differed by " + diff + ", which is greater than the allowed tolerance of " + epsilon);
            VariantContextTestUtils.assertVariantContextsAreEqual(actualVc, expectedVc, Collections.singletonList(infoKey), Collections.emptyList());
        }

        Assert.assertTrue(!expectedVi.hasNext() && !actualVi.hasNext());
    }
}
