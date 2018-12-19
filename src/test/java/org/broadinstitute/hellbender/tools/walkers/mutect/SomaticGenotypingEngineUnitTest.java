package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import static org.testng.Assert.*;

public class SomaticGenotypingEngineUnitTest {

    @Test
    public void testGetGermlineAltAlleleFrequencies() {
        final double defaultAF = 0.001;
        final double nonDefaultAF1 = 0.1;
        final double nonDefaultAF2 = 0.01;
        final Allele Aref = Allele.create("A", true);
        final Allele C = Allele.create("C");
        final Allele G = Allele.create("G");
        final Allele T = Allele.create("T");

        final String source = "SOURCE";
        final int start = 1;
        final int stop = 1;

        //biallelic, vc has the same alt allele
        final List<Allele> altAlleles1 = Arrays.asList(C);
        final VariantContext vc1 = new VariantContextBuilder(source, "1", start, stop, Arrays.asList(Aref, C))
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, new double[] {nonDefaultAF1}).make();
        final double[] af1 = SomaticGenotypingEngine.getGermlineAltAlleleFrequencies(altAlleles1, Optional.of(vc1), defaultAF);
        Assert.assertEquals(af1.length, altAlleles1.size());
        Assert.assertEquals(af1[0], nonDefaultAF1, 0.00001);

        //biallelic, vc has different alt allele
        final List<Allele> altAlleles2 = Arrays.asList(C);
        final VariantContext vc2 = new VariantContextBuilder(source, "1", start, stop, Arrays.asList(Aref, G))
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, new double[] {nonDefaultAF1}).make();
        final double[] af2 = SomaticGenotypingEngine.getGermlineAltAlleleFrequencies(altAlleles2, Optional.of(vc2), defaultAF);
        Assert.assertEquals(af2.length, altAlleles2.size());
        Assert.assertEquals(af2[0], defaultAF, 0.00001);

        //triallelic, same alt alleles
        final List<Allele> altAlleles3 = Arrays.asList(C, G);
        final VariantContext vc3 = new VariantContextBuilder(source, "1", start, stop, Arrays.asList(Aref, C, G))
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, new double[] {nonDefaultAF1, nonDefaultAF2}).make();
        final double[] af3 = SomaticGenotypingEngine.getGermlineAltAlleleFrequencies(altAlleles3, Optional.of(vc3), defaultAF);
        Assert.assertEquals(af3.length, altAlleles3.size());
        Assert.assertEquals(af3[0], nonDefaultAF1, 0.00001);
        Assert.assertEquals(af3[1], nonDefaultAF2, 0.00001);

        //triallelic, same alt alleles in different order
        final List<Allele> altAlleles4 = Arrays.asList(C, G);
        final VariantContext vc4 = new VariantContextBuilder(source, "1", start, stop, Arrays.asList(Aref, G, C))
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, new double[] {nonDefaultAF1, nonDefaultAF2}).make();
        final double[] af4 = SomaticGenotypingEngine.getGermlineAltAlleleFrequencies(altAlleles4, Optional.of(vc4), defaultAF);
        Assert.assertEquals(af4.length, altAlleles4.size());
        Assert.assertEquals(af4[0], nonDefaultAF2, 0.00001);
        Assert.assertEquals(af4[1], nonDefaultAF1, 0.00001);

        //triallelic, only one allele in common
        final List<Allele> altAlleles5 = Arrays.asList(C, G);
        final VariantContext vc5 = new VariantContextBuilder(source, "1", start, stop, Arrays.asList(Aref, C, T))
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, new double[] {nonDefaultAF1, nonDefaultAF2}).make();
        final double[] af5 = SomaticGenotypingEngine.getGermlineAltAlleleFrequencies(altAlleles5, Optional.of(vc5), defaultAF);
        Assert.assertEquals(af5.length, altAlleles5.size());
        Assert.assertEquals(af5[0], nonDefaultAF1, 0.00001);
        Assert.assertEquals(af5[1], defaultAF, 0.00001);
    }

}