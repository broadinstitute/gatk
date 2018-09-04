package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.TextCigarCodec;
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
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public final class CoverageUnitTest extends GATKBaseTest {

    private static final GATKRead makeRead() {
        return ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAllNull() throws Exception {
        final VariantContext vc= null;
        final ReferenceContext referenceContext= null;
        final InfoFieldAnnotation cov = new Coverage();
        final Map<String, Object> annotate = cov.annotate(referenceContext, vc, null); //vc can't be null
    }

    @Test
    public void testDescriptions() throws Exception {
        final InfoFieldAnnotation cov = new Coverage();
        Assert.assertEquals(cov.getDescriptions().size(), 1);
        Assert.assertEquals(cov.getDescriptions().get(0).getID(), VCFConstants.DEPTH_KEY);
    }

    @Test
    public void testLikelihoodsEmpty() throws Exception {
        final List<GATKRead> reads = new ArrayList<>();
        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of("sample1", reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList("sample1"));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(Allele.NO_CALL));
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readsBySample);

        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new Coverage().annotate(referenceContext, vc, likelihoods);
        Assert.assertTrue(annotate.isEmpty());
    }

    private VariantContext makeVC() {
        final GenotypesContext testGC = GenotypesContext.create(2);
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).genotypes(testGC).make();
    }

    @Test
    public void testLikelihoods(){
        final Allele REF = Allele.create("T", true);
        final Allele ALT = Allele.create("A");

        final int refDepth= 5;
        final int altDepth= 3;
        final List<GATKRead>  refReads = IntStream.range(0, refDepth).mapToObj(n -> makeRead()).collect(Collectors.toList());
        final List<GATKRead>  altReads = IntStream.range(0, altDepth).mapToObj(n -> makeRead()).collect(Collectors.toList());
        final ReadLikelihoods<Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods("sample1", refReads, altReads, -100.0, -100.0, REF, ALT);

        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new Coverage().annotate(referenceContext, vc, likelihoods);
        Assert.assertEquals(annotate.size(), 1, "size");
        Assert.assertEquals(annotate.keySet(), Collections.singleton(VCFConstants.DEPTH_KEY), "annots");
        final int m = altDepth + refDepth;
        Assert.assertEquals(annotate.get(VCFConstants.DEPTH_KEY), String.valueOf(m));
    }
}
