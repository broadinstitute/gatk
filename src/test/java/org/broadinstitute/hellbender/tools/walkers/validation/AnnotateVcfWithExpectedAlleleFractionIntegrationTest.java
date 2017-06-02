package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.walkers.validation.AnnotateVcfWithExpectedAlleleFraction;
import org.broadinstitute.hellbender.tools.walkers.validation.MixingFraction;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Created by David Benjamin on 1/31/17.
 */
public class AnnotateVcfWithExpectedAlleleFractionIntegrationTest extends CommandLineProgramTest {

    private static final String VCF_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/mutect/validation/";
    private static final File INPUT_VCF = new File(VCF_DIRECTORY, "dream_4_mixing.vcf");

    // run with made-up mixing fractions and the doctored 2-sample version of DREAM challenge sample 4
    // described in {@link CalculateMixingFractionsIntegrationTest}
    @Test
    public void test() {
        final File table = createTempFile("mixing", ".table");
        final File outputVcf = createTempFile("output", ".vcf");

        final String sample1 = "SAMPLE1";   //as in the input vcf
        final String sample2 = "SAMPLE2";   //as in the input vcf
        final double fraction1 = 0.4;
        final double fraction2 = 0.6;
        MixingFraction.writeMixingFractions(Arrays.asList(new MixingFraction(sample1, fraction1), new MixingFraction(sample2, fraction2)), table);

        final String[] arguments = {
                "-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME, INPUT_VCF.getAbsolutePath(),
                "-" + AnnotateVcfWithExpectedAlleleFraction.MIXING_FRACTIONS_TABLE_NAME, table.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputVcf.getAbsolutePath()
        };

        runCommandLine(arguments);

        final List<VariantContext> input = StreamSupport.stream(new FeatureDataSource<VariantContext>(INPUT_VCF).spliterator(), false)
                .collect(Collectors.toList());
        final List<VariantContext> output = StreamSupport.stream(new FeatureDataSource<VariantContext>(outputVcf).spliterator(), false)
                .collect(Collectors.toList());

        Assert.assertEquals(input.size(), output.size());

        final List<String> inputKeys = input.stream().map(vc -> keyForVariant(vc)).collect(Collectors.toList());
        final List<String> outputKeys = output.stream().map(vc -> keyForVariant(vc)).collect(Collectors.toList());

        Assert.assertEquals(inputKeys, outputKeys);

        final List<Double> alleleFractions = output.stream()
                .map(vc -> vc.getAttributeAsDouble(AnnotateVcfWithExpectedAlleleFraction.EXPECTED_ALLELE_FRACTION_NAME, -1))
                .collect(Collectors.toList());
        // the first few -- 0.2 is sample1 is het, 0.3 is sample 2 is het, 0.5 if both are het
        final List<Double> firstSeveralAlleleFractionsByHand = Arrays.asList(0.2, 0.2, 0.3, 0.2, 0.3, 0.5, 0.3, 0.2, 0.3, 0.3);
        Assert.assertEquals(alleleFractions.subList(0, firstSeveralAlleleFractionsByHand.size()), firstSeveralAlleleFractionsByHand);

        Assert.assertEquals(alleleFractions.get(16), 0.7);  // hom var + het
        Assert.assertEquals(alleleFractions.get(18), 0.8);  // het + hom var
        Assert.assertEquals(alleleFractions.get(26), 0.0);  //both hom ref

    }

    private static String keyForVariant( final VariantContext variant ) {
        return String.format("%s:%d-%d %s", variant.getContig(), variant.getStart(), variant.getEnd(), variant.getAlleles());
    }
}