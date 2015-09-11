package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.Map;

public final class CoverageUnitTest extends BaseTest {

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAllNull() throws Exception {
        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = null;
        final VariantContext vc= null;
        final ReferenceContext referenceContext= null;
        final InfoFieldAnnotation cov = new Coverage();
        final Map<String, Object> annotate = cov.annotate(referenceContext, vc, perReadAlleleLikelihoodMap); //vc can't be null
    }

    @Test
    public void testDescriptions() throws Exception {
        final InfoFieldAnnotation cov = new Coverage();
        Assert.assertEquals(cov.getDescriptions().size(), 1);
        Assert.assertEquals(cov.getDescriptions().get(0).getID(), VCFConstants.DEPTH_KEY);
    }

    @Test
    public void testPerReadAlleleLikelihoodMapEmpty() throws Exception {
        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = Collections.emptyMap();
        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new Coverage().annotate(referenceContext, vc, perReadAlleleLikelihoodMap);
        Assert.assertNull(annotate);
    }

    private VariantContext makeVC() {
        final GenotypesContext testGC = GenotypesContext.create(2);
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).genotypes(testGC).make();
    }

    @Test
    public void testPerReadAlleleLikelihoodMap(){
        final PerReadAlleleLikelihoodMap map= new PerReadAlleleLikelihoodMap();

        final Allele alleleT = Allele.create("T");
        final Allele alleleA = Allele.create("A");
        final double lik= -1.0;

        final int n1A= 3;
        final int n1T= 5;
        for (int i = 0; i < n1A; i++) {
            map.add(ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"), "n1A_" + i), alleleA, lik);
        }
        for (int i = 0; i < n1T; i++) {
            //try to fool it - add 2 alleles for same read
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"), "n1T_" + i);
            map.add(read, alleleA, lik);
            map.add(read, alleleT, lik);
        }

        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = Collections.singletonMap("sample1", map);
        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new Coverage().annotate(referenceContext, vc, perReadAlleleLikelihoodMap);
        Assert.assertEquals(annotate.size(), 1, "size");
        Assert.assertEquals(annotate.keySet(), Collections.singleton(VCFConstants.DEPTH_KEY), "annots");
        final int n= n1A + n1T;
        Assert.assertEquals(annotate.get(VCFConstants.DEPTH_KEY), String.valueOf(n));

    }
}
