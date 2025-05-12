package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.EnvironmentTestUtils;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Iterator;
import java.util.stream.Collectors;

public class NVScoreVariantsIntegrationTest extends CommandLineProgramTest {
    private static final String inputVCF = largeFileTestDir + "VQSR/recalibrated_chr20_start.vcf";
    private static final String inputBAM = largeFileTestDir + "VQSR/g94982_contig_20_start_bamout.bam";
    private static final String reference = b37_reference_20_21;

    private static final double EPSILON_FOR_1D = 0.01;
    private static final double EPSILON_FOR_2D = 0.01;

    @Test(groups = {"python"})
    public void test1DModel() {
        final File tempVcf = createTempFile("test1DModel", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/cnn_1d_chr20_subset_expected.vcf");

        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getAbsolutePath())
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, reference);

        runCommandLine(argsBuilder);

        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_1D_KEY, EPSILON_FOR_1D);
    }

    @Test(groups = {"python"})
    public void test2DModel() {
        final File tempVcf = createTempFile("test2DModel", ".vcf");
        final File expectedVcf = new File(largeFileTestDir + "VQSR/expected/cnn_2d_chr20_subset_expected.vcf");

        final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
        argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getAbsolutePath())
                .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, reference)
                .add(StandardArgumentDefinitions.INPUT_LONG_NAME, inputBAM)
                .add("tensor-type", NVScoreVariants.TensorType.read_tensor.name());

        runCommandLine(argsBuilder);

        assertInfoFieldsAreClose(tempVcf, expectedVcf, GATKVCFConstants.CNN_2D_KEY, EPSILON_FOR_2D);
    }

    @Test(
        expectedExceptions = UserException.NotAvailableInGatkLiteDocker.class,
        singleThreaded = true
    )
    public void testInGatkLiteDocker() {
        EnvironmentTestUtils.checkWithGATKDockerPropertySet(() -> {
            final File tempVcf = createTempFile("test1DModel", ".vcf");

            final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
            argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, inputVCF)
                    .add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempVcf.getAbsolutePath())
                    .add(StandardArgumentDefinitions.REFERENCE_LONG_NAME, reference);

            runCommandLine(argsBuilder);
        });
    }

    private void assertInfoFieldsAreClose(final File actualVcf, final File expectedVcf, final String infoKey, final double epsilon) {
        Iterator<VariantContext> expectedVi = VariantContextTestUtils.streamVcf(expectedVcf).collect(Collectors.toList()).iterator();
        Iterator<VariantContext> actualVi = VariantContextTestUtils.streamVcf(actualVcf).collect(Collectors.toList()).iterator();
        boolean failed = false;
        int totalCount = 0;
        int mismatchCount = 0;
        int annotationMismatchCount = 0;
        while (expectedVi.hasNext() && actualVi.hasNext()) {
            totalCount++;
            VariantContext expectedVc = expectedVi.next();
            VariantContext actualVc = actualVi.next();
            Assert.assertEquals(actualVc.getContig(), expectedVc.getContig(), "Variants from actual and expected VCFs do not match in their location");
            Assert.assertEquals(actualVc.getStart(), expectedVc.getStart(), "Variants from actual and expected VCFs do not match in their location");

            if( expectedVc.hasAttribute(infoKey) != actualVc.hasAttribute(infoKey)) {
                annotationMismatchCount++;
                failed = true;
            }

            double expectedScore = expectedVc.getAttributeAsDouble(infoKey, 0.0); // Different defaults trigger failures on missing scores
            double actualScore = actualVc.getAttributeAsDouble(infoKey, epsilon+1.0);
            double diff = Math.abs(expectedScore-actualScore);
            if ( diff > epsilon) {
                mismatchCount++;
                System.err.println( "scores at " + expectedVc.getContig() + ":" + expectedVc.getStart() + " differed by " + diff
                        + " (expected: " + expectedScore +", actual:" + actualScore + ")," +
                        " which is greater than the allowed tolerance of " + epsilon);
                failed = true;
            }
        }
        if( totalCount == 0 ) {
            failed = true;
        }
        Assert.assertFalse(failed, "Test failed with " + mismatchCount + " significant differences out of " + totalCount +".\n" +
                "There were " + annotationMismatchCount + " sites where the score annotations were present in one file but not the other.");

        Assert.assertTrue(!expectedVi.hasNext() && !actualVi.hasNext());
    }
}
