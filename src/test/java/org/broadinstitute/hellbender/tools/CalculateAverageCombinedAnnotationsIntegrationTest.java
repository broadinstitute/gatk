package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
public class CalculateAverageCombinedAnnotationsIntegrationTest extends CommandLineProgramTest {

    public final String INPUT_VCF = toolsTestDir + "calculate_average_combined_annotations.vcf";
    public final String BAD_INPUT = toolsTestDir + "count_variants.vcf";

    @Test
    public void testExampleVariantWalker() {
        final File outputVcf = createTempFile("calculate_average", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", INPUT_VCF);
        args.add("O", outputVcf);
        args.add(CalculateAverageCombinedAnnotations.ANNOTATION_LIST_SHORT_NAME, GATKVCFConstants.HAPLOTYPES_BEFORE_FILTERING_KEY);
        args.add(CalculateAverageCombinedAnnotations.ANNOTATION_LIST_SHORT_NAME, GATKVCFConstants.HAPLOTYPES_FILTERED_KEY);
        args.add(CalculateAverageCombinedAnnotations.ANNOTATION_LIST_SHORT_NAME, GATKVCFConstants.TREE_SCORE);

        runCommandLine(args);

        VariantContextTestUtils.streamVcf(outputVcf).forEach(vc -> {
            Assert.assertEquals(vc.getAttributeAsDouble("AVERAGE_" + GATKVCFConstants.HAPLOTYPES_BEFORE_FILTERING_KEY, -1), 3.23);
            Assert.assertEquals(vc.getAttributeAsDouble("AVERAGE_" + GATKVCFConstants.HAPLOTYPES_FILTERED_KEY, -1), 1.38);
            Assert.assertEquals(vc.getAttributeAsDouble("AVERAGE_" + GATKVCFConstants.TREE_SCORE, -1), 0.468);
        });
    }

    @Test
    public void testVcfWithWrongAnnotation() {
        final File outputVcf = createTempFile("calculate_average", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("V", BAD_INPUT);
        args.add("O", outputVcf);
        args.add(CalculateAverageCombinedAnnotations.ANNOTATION_LIST_SHORT_NAME, GATKVCFConstants.TREE_SCORE);

        Assert.assertThrows(UserException.class, () -> runCommandLine(args));
    }
}
