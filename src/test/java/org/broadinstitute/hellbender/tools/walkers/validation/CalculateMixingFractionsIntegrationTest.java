package org.broadinstitute.hellbender.tools.walkers.validation;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.walkers.validation.MixingFraction;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Created by David Benjamin on 1/31/17.
 */
public class CalculateMixingFractionsIntegrationTest extends CommandLineProgramTest {
    private static final String VCF_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/mutect/validation/";
    private static final File INPUT_VCF = new File(VCF_DIRECTORY, "dream_4_mixing.vcf");
    private static final String DREAM_BAMS_DIR = largeFileTestDir + "mutect/dream_synthetic_bams/";
    private static final File INPUT_BAM = new File(DREAM_BAMS_DIR, "tumor_4.bam");

    /**
     * The DREAM challenge 4th sample simulates a tumor with subclonal populations of 50% and 35%.  We pretend
     * here that it is actually a mixture of two samples, and we doctored the corresponding such that every variant with allele fraction
     * of roughly 50%/2 = 25% is a het from one sample and every variant with allele fraction roughly 35%/2 is a het from the other
     * sample.  Therefore, the unnormalized mixing fractions should come out to about 0.5 and 0.35, and the normalized fractions
     * should be roughly 0.59 and 0.41.
     *
     * However, when you actually check in IGV, variants cluster into allele fractions of roughly 10% and roughly 25%, which yields
     * normalized mixing fractions of about 30% and 70%, as we in fact obtain.
     */
    @Test
    public void test() {
        final File outputTable = createTempFile("mixing", ".table");

        final String[] arguments = {
                "-" + StandardArgumentDefinitions.VARIANT_SHORT_NAME, INPUT_VCF.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, INPUT_BAM.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputTable.getAbsolutePath()
        };

        runCommandLine(arguments);

        final List<MixingFraction> mixing = MixingFraction.readMixingFractions(outputTable);
        final Map<String, Double> result = mixing.stream().collect(Collectors.toMap(MixingFraction::getSample, MixingFraction::getMixingFraction));
        Assert.assertEquals(result.get("SAMPLE1"), 0.31, 0.02);
        Assert.assertEquals(result.get("SAMPLE2"), 0.69, 0.02);
    }
}