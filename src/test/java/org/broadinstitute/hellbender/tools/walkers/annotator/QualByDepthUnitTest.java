package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.math3.stat.StatUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_QualByDepth;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class QualByDepthUnitTest extends GATKBaseTest {
    private static final Allele REF = Allele.create("A", true);
    private static final Allele ALT = Allele.create("C");

    @DataProvider(name = "UsingAD")
    public Object[][] makeUsingADData() {
        List<Object[]> tests = new ArrayList<>();

        final Allele A = Allele.create("A", true);
        final Allele C = Allele.create("C");
        final Allele G = Allele.create("G");

        final List<Allele> AA = Arrays.asList(A, A);
        final List<Allele> AC = Arrays.asList(A, C);
        final List<Allele> GG = Arrays.asList(G, G);
        final List<Allele> ACG = Arrays.asList(A, C, G);

        final Genotype gAC = new GenotypeBuilder("1", AC).DP(10).AD(new int[]{5,5}).make();
        final Genotype gACNoADNoDP = new GenotypeBuilder("1", AC).AD(new int[]{0,0}).make();
        final Genotype gACNoAD = new GenotypeBuilder("1", AC).DP(10).AD(new int[]{0,0}).make();
        final Genotype gAA = new GenotypeBuilder("2", AA).DP(10).AD(new int[]{10,0}).make();
        final Genotype gACerror = new GenotypeBuilder("3", AC).DP(10).AD(new int[]{9,1}).make();
        final Genotype gACNoError = new GenotypeBuilder("3", AC).DP(17).AD(new int[]{8,2}).make();   //make DP distinctive so that we can test the precedence properly
        final Genotype gGG = new GenotypeBuilder("4", GG).DP(10).AD(new int[]{1,9}).make();

        final double log10PError = -5;
        final double qual = -10.0 * log10PError;
        tests.add(new Object[]{new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make(), qual / 10});   //het
        tests.add(new Object[]{new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).genotypes(Arrays.asList(gACNoADNoDP)).make(), 0.0});   // No AD
        tests.add(new Object[]{new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).genotypes(Arrays.asList(gACNoAD)).make(), qual / 10});   // No AD but present DP, expect it to default to the DP
        tests.add(new Object[]{new VariantContextBuilder("test", "20", 10, 10, ACG).log10PError(log10PError).genotypes(Arrays.asList(gAC, gACNoADNoDP)).make(), qual / 10});   // No AD

        tests.add(new Object[]{new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).genotypes(Arrays.asList(gACerror)).make(), qual / 10});       //het but use 1+X if AD=1,X
        tests.add(new Object[]{new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).genotypes(Arrays.asList(gAA, gAC)).make(), qual / 10});       //het, skip hom
        tests.add(new Object[]{new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).genotypes(Arrays.asList(gAC, gACerror)).make(), qual / 10});         //one het was OK, so use only depth from that one
        tests.add(new Object[]{new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).genotypes(Arrays.asList(gAC, gACNoError)).make(), qual / (10+10)});  //two hets were OK (gACNoError is treated as a regular het)
        tests.add(new Object[]{new VariantContextBuilder("test", "20", 10, 10, ACG).log10PError(log10PError).genotypes(Arrays.asList(gAA, gAC, gACerror, gGG)).make(),  qual / (10 + 10)}); //two hets: gAC and gGG

        final Genotype gACWeirdDP = new GenotypeBuilder("1", AC).DP(70).AD(new int[]{5,5}).make();
        tests.add(new Object[]{new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).genotypes(Arrays.asList(gACWeirdDP)).make(), qual / 10});   //still het, not using DP

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "UsingAD")
    public void testUsingAD(final VariantContext vc, final double expectedQD) {
        final Map<String, Object> annotatedMap = new QualByDepth().annotate(null, vc, null);
        Assert.assertNotNull(annotatedMap, vc.toString());
        if ( annotatedMap.containsKey(GATKVCFConstants.QUAL_BY_DEPTH_KEY) ) {
            final String QD = (String) annotatedMap.get(GATKVCFConstants.QUAL_BY_DEPTH_KEY);
            Assert.assertEquals(Double.valueOf(QD), expectedQD, 0.0001);
            Assert.assertEquals(new QualByDepth().getKeyNames(), Collections.singletonList(GATKVCFConstants.QUAL_BY_DEPTH_KEY));
            Assert.assertEquals(new QualByDepth().getDescriptions(), Collections.singletonList(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.QUAL_BY_DEPTH_KEY)));
        } else {
            Assert.assertTrue(annotatedMap.isEmpty());
        }
    }

    @Test
    public void testUsingDP() {
        final List<Allele> AC = Arrays.asList(REF, ALT);
        final int depth = 20;
        final Genotype gAC = new GenotypeBuilder("1", AC).DP(depth).make();

        final double log10PError = -5;
        final double qual = -10.0 * log10PError;

        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();

        final Map<String, Object> annotatedMap = new QualByDepth().annotate(null, vc, null);
        Assert.assertNotNull(annotatedMap, vc.toString());
        final String QD = (String)annotatedMap.get(GATKVCFConstants.QUAL_BY_DEPTH_KEY);

        final double expectedQD = qual/depth;
        Assert.assertEquals(Double.valueOf(QD), expectedQD, 0.0001);
    }

    @Test
    public void testUsingReads(){
        final List<Allele> ALLELES = Arrays.asList(REF, ALT);
        final int depth = 20;
        final String sample1 = "sample1";
        final int dpDepth = 30; //Note: using a different value on purpose so that we can check that reads are preferred over DP
        final Genotype gAC = new GenotypeBuilder(sample1, ALLELES).DP(dpDepth).make();

        final double log10PError = -5;
        final double qual = -10.0 * log10PError;

        final List<GATKRead> reads = IntStream.range(0, depth)
                .mapToObj(n -> ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"))).collect(Collectors.toList());

        final ReadLikelihoods<Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods(sample1, reads, -100.0, REF, ALT);

        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, ALLELES).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();
        final Map<String, Object> annotatedMap = new QualByDepth().annotate(null, vc, likelihoods);
        Assert.assertNotNull(annotatedMap, vc.toString());
        final String QD = (String)annotatedMap.get(GATKVCFConstants.QUAL_BY_DEPTH_KEY);

        final double expectedQD = qual/depth;
        Assert.assertEquals(Double.valueOf(QD), expectedQD, 0.0001);

        //Now we test that when AD is present, it trumps everything
        final Genotype gAC_withAD = new GenotypeBuilder("1", ALLELES).DP(dpDepth).AD(new int[]{5,5}).make();
        final VariantContext vc_withAD = new VariantContextBuilder("test", "20", 10, 10, ALLELES).log10PError(log10PError).genotypes(Arrays.asList(gAC_withAD)).make();
        final Map<String, Object> annotatedMap_withAD = new QualByDepth().annotate(null, vc_withAD, likelihoods);
        final String QD_withAD = (String)annotatedMap_withAD.get(GATKVCFConstants.QUAL_BY_DEPTH_KEY);
        final double expectedQD_withAD = qual/(5+5);//two AD fields
        Assert.assertEquals(Double.valueOf(QD_withAD), expectedQD_withAD, 0.0001);

    }

    @DataProvider(name = "nullQD")
    public Object[][] nullQD() {
        List<Object[]> tests = new ArrayList<>();

        final Allele A = Allele.create("A", true);
        final Allele C = Allele.create("C");

        final List<Allele> AA = Arrays.asList(A, A);
        final List<Allele> AC = Arrays.asList(A, C);

        final Genotype gAC = new GenotypeBuilder("1", AC).DP(10).AD(new int[]{5,5}).make();
        final Genotype gAA = new GenotypeBuilder("2", AA).DP(10).AD(new int[]{10,0}).make();

        final double log10PError = -5;
        tests.add(new Object[]{new VariantContextBuilder("test", "20", 10, 10, AC).genotypes(Arrays.asList(gAC)).make()});   //no error
        tests.add(new Object[]{new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).make()}); //no genotypes
        tests.add(new Object[]{new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).genotypes(Arrays.asList(gAA)).make()}); //no hets

        return tests.toArray(new Object[][]{});
    }
    @Test(dataProvider = "nullQD")
    public void testEmpty(final VariantContext vc){
        final Map<String, Object> annotatedMap = new QualByDepth().annotate(null, vc, null);
        Assert.assertTrue(annotatedMap.isEmpty());
    }

    @Test
    public void testVeryHighQD(){
        final Allele A = Allele.create("A", true);
        final Allele C = Allele.create("C");

        final List<Allele> AC = Arrays.asList(A, C);
        final Genotype gAC = new GenotypeBuilder("1", AC).DP(10).AD(new int[]{5, 5}).make();

        final double lowError = -QualByDepth.MAX_QD_BEFORE_FIXING - 10;

        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(lowError).genotypes(Arrays.asList(gAC)).make();

        final double[] qds = new double[100000];
        for (int i = 0; i < qds.length; i++) {
            final Map<String, Object> annotatedMap = new QualByDepth().annotate(null, vc, null);
            final String QD = (String)annotatedMap.get(GATKVCFConstants.QUAL_BY_DEPTH_KEY);
            final double qdVal = Double.valueOf(QD);
            qds[i] = qdVal;
        }
        //test that on average it'll get settled at the expected mean
        Assert.assertEquals(StatUtils.mean(qds), QualByDepth.IDEAL_HIGH_QD, 0.02);
    }

    @Test
    public void testAnnotate_AS() throws Exception {

        final Allele A = Allele.create("A", true);
        final Allele C = Allele.create("C");
        final Allele G = Allele.create("G");

        final List<Allele> AC = Arrays.asList(A, C);
        final int depth = 20;
        final String sample1 = "sample1";
        final int dpDepth = 30; //Note: using a different value on purpose so that we can check that reads are preferred over DP
        final Genotype gAC = new GenotypeBuilder(sample1, AC).DP(dpDepth).make();

        final double log10PError = -5;

        final List<GATKRead> reads = IntStream.range(0, depth)
                .mapToObj(n -> ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"))).collect(Collectors.toList());

        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of(sample1, reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList(sample1));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(A, C, G));
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readsBySample);

        // modify likelihoods in-place
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(0);

        for (int n = 0; n < depth; n++) {
            matrix.set(0, n, -1.0);
            matrix.set(1, n, -100.0);
            matrix.set(2, n, -1000.0);
        }

        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();
        final Map<String, Object> annotatedMap = new AS_QualByDepth().annotate(null, vc, likelihoods);
        Assert.assertTrue(annotatedMap.isEmpty());

        final Map<String, Object> annotatedMapRaw = new AS_QualByDepth().annotateRawData(null, vc, likelihoods);
        Assert.assertTrue(annotatedMapRaw.isEmpty());
    }

    @Test
    public void testDescriptions_AS() throws Exception {
        final AS_QualByDepth cov = new AS_QualByDepth();
        Assert.assertEquals(cov.getRawKeyName(), GATKVCFConstants.AS_QUAL_KEY);
        Assert.assertEquals(cov.getDescriptions().size(), 1);
        Assert.assertEquals(cov.getDescriptions().get(0).getID(), GATKVCFConstants.AS_QUAL_BY_DEPTH_KEY);
    }
}
