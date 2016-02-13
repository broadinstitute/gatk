package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.Map;

public final class MappingQualityZeroUnitTest {

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAllNull() throws Exception {
        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = null;
        final VariantContext vc= null;
        final ReferenceContext referenceContext= null;
        final InfoFieldAnnotation cov = new MappingQualityZero();
        final Map<String, Object> annotate = cov.annotate(referenceContext, vc, perReadAlleleLikelihoodMap);//vc can't be null
    }

    @Test
    public void testNullStratifiedPerReadAlleleLikelihoodMap() throws Exception {
        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = null;
        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final InfoFieldAnnotation cov = new MappingQualityZero();
        final Map<String, Object> annotate = cov.annotate(referenceContext, vc, perReadAlleleLikelihoodMap);
        Assert.assertNull(annotate);

        Assert.assertEquals(cov.getDescriptions().size(), 1);
        Assert.assertEquals(cov.getDescriptions().get(0).getID(), VCFConstants.MAPPING_QUALITY_ZERO_KEY);
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
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"), "n1A_" + i);
            read.setMappingQuality(10);
            map.add(read, alleleA, lik);
        }
        for (int i = 0; i < n1T; i++) {
            //try to fool it - add 2 alleles for same read
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"), "n1T_" + i);
            read.setMappingQuality(0);
            map.add(read, alleleA, lik);
            map.add(read, alleleT, lik);
        }

        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = Collections.singletonMap("sample1", map);
        final Map<String, AlignmentContext> stratifiedContexts= null;
        final VariantContext vc = makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new MappingQualityZero().annotate(referenceContext, vc, perReadAlleleLikelihoodMap);
        Assert.assertEquals(annotate.size(), 1, "size");
        Assert.assertEquals(annotate.keySet(), Collections.singleton(VCFConstants.MAPPING_QUALITY_ZERO_KEY), "annots");
        final int n= n1T; //only those are MQ0
        Assert.assertEquals(annotate.get(VCFConstants.MAPPING_QUALITY_ZERO_KEY), String.valueOf(n));

    }


    @Test
    public void testPerReadAlleleLikelihoodMapEmpty() throws Exception {
        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = Collections.emptyMap();
        final Map<String, AlignmentContext> stratifiedContexts= null;
        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new MappingQualityZero().annotate(referenceContext, vc, perReadAlleleLikelihoodMap);

        final int n= 0; //strangely,  MappingQualityZero returns 0 if perReadAlleleLikelihoodMap is empty
        Assert.assertEquals(annotate.get(VCFConstants.MAPPING_QUALITY_ZERO_KEY), String.valueOf(n));
    }
}
