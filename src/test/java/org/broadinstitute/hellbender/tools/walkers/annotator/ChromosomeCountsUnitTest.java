package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Map;

public final class ChromosomeCountsUnitTest {
    private static final Allele REF = Allele.create("A", true);
    private static final Allele ALT = Allele.create("T");

    private static final String SAMPLE_1 = "NA1";
    private static final String SAMPLE_2 = "NA2";

    private VariantContext makeVC() {
        final double[] genotypeLikelihoods1 = {30,0,190};
        final GenotypesContext testGC = GenotypesContext.create(2);
        // SAMPLE_1 -> A/T with GQ 30
        testGC.add(new GenotypeBuilder(SAMPLE_1).alleles(Arrays.asList(REF, ALT)).PL(genotypeLikelihoods1).GQ(30).make());
        // SAMPLE_2 -> A/T with GQ 40
        testGC.add(new GenotypeBuilder(SAMPLE_2).alleles(Arrays.asList(REF, ALT)).PL(genotypeLikelihoods1).GQ(40).make());

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(REF, ALT)).chr("1").start(15L).stop(15L).genotypes(testGC).make();
    }

    private VariantContext makeEmptyVC() {
        return new VariantContextBuilder().alleles(Arrays.asList(REF, ALT)).chr("1").start(15L).stop(15L).make();
    }

    @Test
    public void testVC() throws Exception {
        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final InfoFieldAnnotation ann = new ChromosomeCounts();
        final Map<String, Object> annotate = ann.annotate(referenceContext, vc, null);

        //two hets
        Assert.assertEquals(annotate.get(VCFConstants.ALLELE_NUMBER_KEY), 4);
        Assert.assertEquals(annotate.get(VCFConstants.ALLELE_COUNT_KEY), 2);
        Assert.assertEquals(annotate.get(VCFConstants.ALLELE_FREQUENCY_KEY), 0.5);

        Assert.assertEquals(ann.getDescriptions(), Arrays.asList(ChromosomeCounts.descriptions));
        Assert.assertEquals(ann.getKeyNames(), Arrays.asList(ChromosomeCounts.keyNames));
    }
    
    @Test
    public void testEmptyVC() throws Exception {
        final VariantContext vc = makeEmptyVC();
        final ReferenceContext referenceContext = null;
        final InfoFieldAnnotation cov = new ChromosomeCounts();
        final Map<String, Object> annotate = cov.annotate(referenceContext, vc, null);

        //no genotypes
        Assert.assertTrue(annotate.isEmpty());
    }
}
