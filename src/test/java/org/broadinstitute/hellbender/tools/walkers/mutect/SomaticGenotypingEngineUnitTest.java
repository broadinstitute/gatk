package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import static org.testng.Assert.*;

public class SomaticGenotypingEngineUnitTest {

    @Test
    public void testGetGermlineAltAlleleFrequencies() {
        final double defaultAF = 0.001;
        final double nonDefaultAF1 = 0.1;
        final double nonDefaultAF2 = 0.01;

        final Allele Aref = Allele.create("A", true);
        final Allele CTref = Allele.create("CT", true);
        final Allele Cref = Allele.create("C", true);

        final Allele C = Allele.create("C");
        final Allele G = Allele.create("G");
        final Allele T = Allele.create("T");
        final Allele CT = Allele.create("CT");
        final Allele CTT = Allele.create("CTT");
        final Allele CTTT = Allele.create("CTTT");

        final String source = "SOURCE";
        final int start = 1;
        final int stop = 1;

        //biallelic, vc has the same alt allele
        final VariantContext vc1 = new VariantContextBuilder(source, "1", start, stop, Arrays.asList(Aref, C))
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, new double[] {nonDefaultAF1}).make();
        final List<VariantContext> germlineVC1 = Arrays.asList(vc1);
        final VariantContext emitVC1 = new VariantContextBuilder(source, "1", start, stop, Arrays.asList(Aref, C)).make();
        final Map<Allele, Double> af1 = SomaticGenotypingEngine.getPopulationAFAnnotation(germlineVC1, emitVC1, defaultAF);
        Assert.assertEquals(af1.size(), 1);
        Assert.assertEquals(af1.get(C), nonDefaultAF1, 0.00001);

        //biallelic, vc has different alt allele
        final VariantContext vc2 = new VariantContextBuilder(source, "1", start, stop, Arrays.asList(Aref, G))
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, new double[] {nonDefaultAF1}).make();
        final List<VariantContext> germlineVC2 = Arrays.asList(vc2);
        final VariantContext emitVC2 = new VariantContextBuilder(source, "1", start, stop, Arrays.asList(Aref, C)).make();
        final Map<Allele, Double> af2 = SomaticGenotypingEngine.getPopulationAFAnnotation(germlineVC2, emitVC2, defaultAF);
        Assert.assertEquals(af2.size(), 1);
        Assert.assertEquals(af2.get(C), defaultAF, 0.00001);

        //triallelic, same alt alleles
        final VariantContext vc3 = new VariantContextBuilder(source, "1", start, stop, Arrays.asList(Aref, C, G))
            .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, new double[] {nonDefaultAF1, nonDefaultAF2}).make();
        final List<VariantContext> germlineVC3 = Arrays.asList(vc3);
        final VariantContext emitVC3 = new VariantContextBuilder(source, "1", start, stop, Arrays.asList(Aref, C, G)).make();
        final Map<Allele, Double> af3 = SomaticGenotypingEngine.getPopulationAFAnnotation(germlineVC3, emitVC3, defaultAF);
        Assert.assertEquals(af3.size(), 2);
        Assert.assertEquals(af3.get(C), nonDefaultAF1, 0.00001);
        Assert.assertEquals(af3.get(G), nonDefaultAF2, 0.00001);

        //triallelic, same alt alleles in different order
        final VariantContext vc4 = new VariantContextBuilder(source, "1", start, stop, Arrays.asList(Aref, C, G))
            .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, new double[] {nonDefaultAF1, nonDefaultAF2}).make();
        final List<VariantContext> germlineVC4 = Arrays.asList(vc4);
        final VariantContext emitVC4 = new VariantContextBuilder(source, "1", start, stop, Arrays.asList(Aref, G, C)).make();
        final Map<Allele, Double> af4 = SomaticGenotypingEngine.getPopulationAFAnnotation(germlineVC4, emitVC4, defaultAF);
        Assert.assertEquals(af4.size(), 2);
        Assert.assertEquals(af4.get(C), nonDefaultAF1, 0.00001);
        Assert.assertEquals(af4.get(G), nonDefaultAF2, 0.00001);

        //triallelic, only one allele in common
        final VariantContext vc5 = new VariantContextBuilder(source, "1", start, stop, Arrays.asList(Aref, C, T))
            .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, new double[] {nonDefaultAF1, nonDefaultAF2}).make();
        final List<VariantContext> germlineVC5 = Arrays.asList(vc5);
        final VariantContext emitVC5 = new VariantContextBuilder(source, "1", start, stop, Arrays.asList(Aref, C, G)).make();
        final Map<Allele, Double> af5 = SomaticGenotypingEngine.getPopulationAFAnnotation(germlineVC5, emitVC5, defaultAF);
        Assert.assertEquals(af5.size(), 2);
        Assert.assertEquals(af5.get(C), nonDefaultAF1, 0.00001);
        Assert.assertEquals(af5.get(G), defaultAF, 0.00001);

        //triallelic reference with mixed DEL/INS/SNP, e.g. CT - C/CTT
        final VariantContext vc6 = new VariantContextBuilder(source, "1", start, stop+1, Arrays.asList(CTref, C, CTT))
            .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, new double[] {nonDefaultAF1, nonDefaultAF2}).make();
        final List<VariantContext> germlineVC6 = Arrays.asList(vc6);
        final VariantContext emitVC6 = new VariantContextBuilder(source, "1", start, stop, Arrays.asList(Cref, CT, CTT)).make();
        final Map<Allele, Double> af6 = SomaticGenotypingEngine.getPopulationAFAnnotation(germlineVC6, emitVC6, defaultAF);
        Assert.assertEquals(af6.size(), 2);
        Assert.assertEquals(af6.get(CT), nonDefaultAF2, 0.00001); //Cref/CT matches CTref/CTT
        Assert.assertEquals(af6.get(CTT), defaultAF, 0.00001);    //No match

        //triallelic to-be-emitted call with mixed DEL/INS/SNP, e.g. CT - C/CTT
        final VariantContext vc7 = new VariantContextBuilder(source, "1", start, stop, Arrays.asList(Cref, CT, CTT))
            .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, new double[] {nonDefaultAF1, nonDefaultAF2}).make();
        final List<VariantContext> germlineVC7 = Arrays.asList(vc7);
        final VariantContext emitVC7 = new VariantContextBuilder(source, "1", start, stop+1, Arrays.asList(CTref, C, CTT)).make();
        final Map<Allele, Double> af7 = SomaticGenotypingEngine.getPopulationAFAnnotation(germlineVC7, emitVC7, defaultAF);
        Assert.assertEquals(af7.size(), 2);
        Assert.assertEquals(af7.get(C), defaultAF, 0.00001);          // no match
        Assert.assertEquals(af7.get(CTT), nonDefaultAF1, 0.00001);    // CTref/CTT matches Cref/CT
    }

}
