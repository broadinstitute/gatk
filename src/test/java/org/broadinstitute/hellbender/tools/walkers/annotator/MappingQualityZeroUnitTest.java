package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public final class MappingQualityZeroUnitTest {

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAllNull() throws Exception {
        final VariantContext vc= null;
        final ReferenceContext referenceContext= null;
        final InfoFieldAnnotation cov = new MappingQualityZero();
        final Map<String, Object> annotate = cov.annotate(referenceContext, vc, null);//vc can't be null
    }

    @Test
    public void testNullLikelihoods() throws Exception {
        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final InfoFieldAnnotation cov = new MappingQualityZero();
        final Map<String, Object> annotate = cov.annotate(referenceContext, vc, null);
        Assert.assertTrue(annotate.isEmpty());

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

    private GATKRead makeRead(final int mappingQual) {
        return ArtificialAnnotationUtils.makeRead(25, mappingQual);
    }

    @Test
    public void testLikelihoods(){

        final Allele REF = Allele.create("T", true);
        final Allele ALT = Allele.create("A");

        final int refDepth= 5;
        final int altDepth= 3;
        final int refMQ = 10;
        final int altMQ = 0;
        final List<GATKRead> refReads = IntStream.range(0, refDepth).mapToObj(i -> makeRead(refMQ)).collect(Collectors.toList());
        final List<GATKRead> altReads = IntStream.range(0, altDepth).mapToObj(i -> makeRead(altMQ)).collect(Collectors.toList());
        final ReadLikelihoods<Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods("sample1", refReads, altReads, -1.0, -1.0, REF, ALT);

        final VariantContext vc = makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new MappingQualityZero().annotate(referenceContext, vc, likelihoods);
        Assert.assertEquals(annotate.size(), 1, "size");
        Assert.assertEquals(annotate.keySet(), Collections.singleton(VCFConstants.MAPPING_QUALITY_ZERO_KEY), "annots");
        Assert.assertEquals(annotate.get(VCFConstants.MAPPING_QUALITY_ZERO_KEY), String.valueOf(altDepth));

    }

    @Test
    public void testEmptyLikelihoods() throws Exception {
        final List<GATKRead> reads = Collections.emptyList();
        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of("sample1", reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList("sample1"));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(Allele.create("A")));
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readsBySample);

        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new MappingQualityZero().annotate(referenceContext, vc, likelihoods);

        final int n= 0; //strangely,  MappingQualityZero returns 0 if likelihoods is empty
        Assert.assertEquals(annotate.get(VCFConstants.MAPPING_QUALITY_ZERO_KEY), String.valueOf(n));
    }
}
