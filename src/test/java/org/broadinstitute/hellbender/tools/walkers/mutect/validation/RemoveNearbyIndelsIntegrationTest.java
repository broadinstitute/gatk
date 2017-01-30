package org.broadinstitute.hellbender.tools.walkers.mutect.validation;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.Main;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.exome.segmentation.PerformAlleleFractionSegmentation;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import static org.testng.Assert.*;

/**
 * Created by davidben on 1/31/17.
 */
public class RemoveNearbyIndelsIntegrationTest extends CommandLineProgramTest {

    private static final String TOOLS_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/mutect/validation/";
    private static final File INPUT_VCF = new File(TOOLS_TEST_DIRECTORY, "nearby_indels.vcf");
    /**
     * nearby_indels.vcf looks like this:
     * #CHROM  POS     ID      REF     ALT
     * 20      48      .       C       GA
     * 20      53      .       T       AT
     * 20      59      .       G       C
     * 20      69      .       C       AG
     * 20      71      .       G       T
     * 20      77      .       A       G
     * 20      89      .       G       A
     * 20      95      .       T       CC
     * 20      100     .       C       A
     *
     * Reoving nearby indels with a threshold of 20 bp should yield:
     * #CHROM  POS     ID      REF     ALT
     * 20      59      .       G       C
     * 20      71      .       G       T
     * 20      77      .       A       G
     * 20      89      .       G       A
     * 20      95      .       T       CC
     * 20      100     .       C       A
     */
    @Test
    public void test() {
        final File inputVcf = INPUT_VCF;

        final File outputVcf = createTempFile("filtered", ".vcf");
        final int minIndelDistance = 20;
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME, inputVcf.getAbsolutePath(),
                "-" + RemoveNearbyIndels.MIN_INDEL_SPACING_NAME, Integer.toString(minIndelDistance),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputVcf.getAbsolutePath()
        };

        runCommandLine(arguments);

        final List<VariantContext> output = StreamSupport.stream(new FeatureDataSource<VariantContext>(outputVcf).spliterator(), false)
                .collect(Collectors.toList());

        final List<Integer> outputPositions = output.stream().map(VariantContext::getStart).collect(Collectors.toList());
        final List<Integer> expectedOutputPositions = Arrays.asList(59, 71, 77, 89, 95, 100);

        Assert.assertEquals(outputPositions, expectedOutputPositions);
    }
}