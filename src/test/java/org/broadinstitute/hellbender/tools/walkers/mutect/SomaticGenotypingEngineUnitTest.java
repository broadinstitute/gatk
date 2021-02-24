package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

public class SomaticGenotypingEngineUnitTest {
    private static final double DEFAULT_AF = 1e-7;

    private static final Allele ATref = Allele.create("AT",true);
    private static final Allele GT = Allele.create("GT",false);

    @DataProvider
    Object[][] germlineAltAlleleFrequenciesData() {
        // {called alleles, germline alleles, germline alt AFs, expected}
        return new Object[][] {
                //biallelic, same alt allele
                { Arrays.asList(Allele.REF_A, Allele.ALT_C), Arrays.asList(Allele.REF_A, Allele.ALT_C), new double[] {0.1}, new double[] {0.1}},

                //biallelic, different alt allele
                { Arrays.asList(Allele.REF_A, Allele.ALT_C), Arrays.asList(Allele.REF_A, Allele.ALT_G), new double[] {0.1}, new double[] {DEFAULT_AF}},

                //triallelic, same alt alleles in same order
                { Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_G), Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_G), new double[] {0.1, 0.2}, new double[] {0.1, 0.2}},

                //triallelic, same alt alleles in different order
                { Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_G), Arrays.asList(Allele.REF_A, Allele.ALT_G, Allele.ALT_C), new double[] {0.1, 0.2}, new double[] {0.2, 0.1}},

                //triallelic, only one alt in common
                { Arrays.asList(Allele.REF_A, Allele.ALT_C, Allele.ALT_G), Arrays.asList(Allele.REF_A, Allele.ALT_G, Allele.ALT_T), new double[] {0.1, 0.2}, new double[] {DEFAULT_AF, 0.1}},

                //biallelic, same alt allele in different representations
                { Arrays.asList(Allele.REF_A, Allele.ALT_G), Arrays.asList(ATref, GT), new double[] {0.1}, new double[] {0.1}},
                { Arrays.asList(ATref, GT), Arrays.asList(Allele.REF_A, Allele.ALT_G), new double[] {0.1}, new double[] {0.1}},
        };
    }

    @Test(dataProvider = "germlineAltAlleleFrequenciesData")
    public void testGetGermlineAltAlleleFrequencies(final List<Allele> calledAlleles, final List<Allele> germlineAlleles, final double[] germlineAltAFs, final double[] expected) {
        final int start = 1;
        final int length = germlineAlleles.stream().mapToInt(Allele::length).max().getAsInt();
        final int stop = start + length - 1;
        final VariantContext vc1 = new VariantContextBuilder("SOURCE", "1", start, stop, germlineAlleles)
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, germlineAltAFs).make();
        final double[] result = SomaticGenotypingEngine.getGermlineAltAlleleFrequencies(calledAlleles, Optional.of(vc1), DEFAULT_AF);
        Assert.assertEquals(result, expected, 1.0e-10);
    }
}