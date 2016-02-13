package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;

public final class RMSMappingQualityUnitTest {

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAllNull() throws Exception {
        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = null;
        final VariantContext vc= null;
        final ReferenceContext referenceContext= null;
        final InfoFieldAnnotation cov = new RMSMappingQuality();
        final Map<String, Object> annotate = cov.annotate(referenceContext, vc, perReadAlleleLikelihoodMap); //vc can't be null
    }

    @Test
    public void testDescriptions() throws Exception {
        final InfoFieldAnnotation cov = new RMSMappingQuality();
        Assert.assertEquals(cov.getDescriptions().size(), 1);
        Assert.assertEquals(cov.getDescriptions().get(0).getID(), VCFConstants.RMS_MAPPING_QUALITY_KEY);
    }

    @Test
    public void testNullStratifiedPerReadAlleleLikelihoodMap() throws Exception {
        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = null;
        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final InfoFieldAnnotation cov = new RMSMappingQuality();
        final Map<String, Object> annotate = cov.annotate(referenceContext, vc, perReadAlleleLikelihoodMap);
        Assert.assertNull(annotate);

        Assert.assertEquals(cov.getDescriptions().size(), 1);
        Assert.assertEquals(cov.getDescriptions().get(0).getID(), VCFConstants.RMS_MAPPING_QUALITY_KEY);
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

        final Allele alleleA = Allele.create("A");
        final double lik= -1.0;  //ignored

        final int[] MQs = {1,2,3,4,5,6,7,8,9,10, QualityUtils.MAPPING_QUALITY_UNAVAILABLE};
        final List<Integer> MQsList = Arrays.asList(ArrayUtils.toObject(MQs));

        //MQ 255 are excluded from the calculations, we test it here.
        final List<Integer> MQsListOK = new ArrayList<>(MQsList);
        //NOTE: if we just call remove(i), Java thinks i is an index.
        //A workaround for this overloading bogosity to to call removeAll and pass a collection
        //(casting i to (Object) would work too but it's more error prone)
        MQsListOK.removeAll(Collections.singleton(QualityUtils.MAPPING_QUALITY_UNAVAILABLE));

        final int n1A= MQs.length;
        for (int i = 0; i < n1A; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"));
            read.setMappingQuality(MQs[i]);
            map.add(read, alleleA, lik);
        }

        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = Collections.singletonMap("sample1", map);
        final VariantContext vc = makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new RMSMappingQuality().annotate(referenceContext, vc, perReadAlleleLikelihoodMap);
        Assert.assertEquals(annotate.size(), 1, "size");
        Assert.assertEquals(annotate.keySet(), Collections.singleton(VCFConstants.RMS_MAPPING_QUALITY_KEY), "annots");
        final double rms= MathUtils.rms(MQsListOK); //only those are MQ0
        Assert.assertEquals(annotate.get(VCFConstants.RMS_MAPPING_QUALITY_KEY), String.format("%.2f", rms));
    }

    @Test
    public void testPerReadAlleleLikelihoodMapEmpty() throws Exception {
        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = Collections.emptyMap();
        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new RMSMappingQuality().annotate(referenceContext, vc, perReadAlleleLikelihoodMap);
        Assert.assertNull(annotate);
    }
}
