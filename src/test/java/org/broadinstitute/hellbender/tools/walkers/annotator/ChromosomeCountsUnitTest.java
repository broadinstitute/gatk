package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Map;

public final class ChromosomeCountsUnitTest {

    private final String sample1 = "NA1";
    private final String sample2 = "NA2";

    private VariantContext makeVC() {
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");
        final double[] genotypeLikelihoods1 = {30,0,190};
        final GenotypesContext testGC = GenotypesContext.create(2);
        // sample1 -> A/T with GQ 30
        testGC.add(new GenotypeBuilder(sample1).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(30).make());
        // sample2 -> A/T with GQ 40
        testGC.add(new GenotypeBuilder(sample2).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(40).make());

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).genotypes(testGC).make();
    }

    private VariantContext makeEmptyVC() {
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).make();
    }

    @Test
    public void testVC() throws Exception {
        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = null;
        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final InfoFieldAnnotation ann = new ChromosomeCounts();
        final Map<String, Object> annotate = ann.annotate(referenceContext, vc, perReadAlleleLikelihoodMap);

        //two hets
        Assert.assertEquals(annotate.get(VCFConstants.ALLELE_NUMBER_KEY), 4);
        Assert.assertEquals(annotate.get(VCFConstants.ALLELE_COUNT_KEY), 2);
        Assert.assertEquals(annotate.get(VCFConstants.ALLELE_FREQUENCY_KEY), 0.5);

        Assert.assertEquals(ann.getDescriptions(), Arrays.asList(ChromosomeCounts.descriptions));
        Assert.assertEquals(ann.getKeyNames(), Arrays.asList(ChromosomeCounts.keyNames));
    }


    @Test
    public void testEmptyVC() throws Exception {
        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = null;
        final VariantContext vc = makeEmptyVC();
        final ReferenceContext referenceContext = null;
        final InfoFieldAnnotation cov = new ChromosomeCounts();
        final Map<String, Object> annotate = cov.annotate(referenceContext, vc, perReadAlleleLikelihoodMap);

        //no genotypes
        Assert.assertNull(annotate);
    }
}
