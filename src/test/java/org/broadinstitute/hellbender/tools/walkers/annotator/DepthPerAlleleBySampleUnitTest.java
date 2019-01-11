package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public final class DepthPerAlleleBySampleUnitTest extends GATKBaseTest {

    private static final Allele REF = Allele.create("A", true);
    private static final Allele ALT = Allele.create("C");
    private static final List<Allele> ALLELES = Arrays.asList(REF, ALT);
    private static final String SAMPLE = "sample1";

    private GATKRead makeRead() {
        return ArtificialAnnotationUtils.makeRead(30, 50);
    }

    @Test
    public void testDescription(){
        Assert.assertEquals(new DepthPerAlleleBySample().getKeyNames(), Collections.singletonList(VCFConstants.GENOTYPE_ALLELE_DEPTHS));
        Assert.assertEquals(new DepthPerAlleleBySample().getDescriptions(), Collections.singletonList(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS)));
    }

    @DataProvider(name = "ReadDepthData")
    public Object[][] readDepthData() {
        return new Object[][]{
                {20, 17},
                {0, 0}
        };
    }

    @Test(dataProvider = "ReadDepthData")
    public void testUsingReads(final int refDepth, final int altDepth){
        final int[] expectedAD = {refDepth, altDepth};

        final int dpDepth = 30; //Note: using a different value on purpose so that we can check that reads are preferred over DP
        final Genotype gAC = new GenotypeBuilder(SAMPLE, ALLELES).DP(dpDepth).make();

        final double log10PError = -5;

        final List<GATKRead> refReads = IntStream.range(0, refDepth).mapToObj(i -> makeRead()).collect(Collectors.toList());
        final List<GATKRead> altReads = IntStream.range(0, altDepth).mapToObj(i -> makeRead()).collect(Collectors.toList());
        final ReadLikelihoods<Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods(SAMPLE, refReads, altReads, -100.0, -100.0, REF, ALT);

        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, ALLELES).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();

        final GenotypeBuilder gb = new GenotypeBuilder(gAC);
        new DepthPerAlleleBySample().annotate(null, vc, gAC, gb, likelihoods);
        final int[] ad = gb.make().getAD();
        Assert.assertEquals(ad, expectedAD);

        //now test a no-op
        final GenotypeBuilder gb1 = new GenotypeBuilder(gAC);
        new DepthPerAlleleBySample().annotate(null, vc, null, gb1, likelihoods);  //null genotype
        Assert.assertFalse(gb1.make().hasAD());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBlowUp(){
        final int dpDepth = 30; //Note: using a different value on purpose so that we can check that reads are preferred over DP
        final Genotype gAC = new GenotypeBuilder(SAMPLE, ALLELES).DP(dpDepth).make();

        final double log10PError = -5;

        final List<GATKRead> reads = Arrays.asList(ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M")));
        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of(SAMPLE, reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList(SAMPLE));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(REF));
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readsBySample);
        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, ALLELES).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();

        final GenotypeBuilder gb = new GenotypeBuilder(gAC);
        //this blows up because there's no C allele in the likelihoods
        new DepthPerAlleleBySample().annotate(null, vc, gAC, gb, likelihoods);
    }

    @Test
    public void testEmptyLikelihoodsFallback(){
        final int dpDepth = 30;
        final Genotype gAC = new GenotypeBuilder(SAMPLE, ALLELES).DP(dpDepth).make();

        final double log10PError = -5;

        Map<String, List<GATKRead>> emptyMap = new HashMap<>();
        emptyMap.put(SAMPLE, Collections.emptyList());
        final ReadLikelihoods<Allele> likelihoods = new UnfilledReadsLikelihoods<Allele>(new IndexedSampleList(SAMPLE), new IndexedAlleleList<>(ALLELES), emptyMap);

        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, ALLELES).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();

        final GenotypeBuilder gb = new GenotypeBuilder(gAC);
        new DepthPerAlleleBySample().annotate(null, vc, gAC, gb, likelihoods);
        final int[] ad = gb.make().getAD();
        Assert.assertEquals(ad, null);

    }

}
