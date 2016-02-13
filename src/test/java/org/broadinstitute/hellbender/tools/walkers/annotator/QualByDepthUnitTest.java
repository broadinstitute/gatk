package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.math3.stat.StatUtils;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class QualByDepthUnitTest extends BaseTest {

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
        final Genotype gAA = new GenotypeBuilder("2", AA).DP(10).AD(new int[]{10,0}).make();
        final Genotype gACerror = new GenotypeBuilder("3", AC).DP(10).AD(new int[]{9,1}).make();
        final Genotype gACNoError = new GenotypeBuilder("3", AC).DP(17).AD(new int[]{8,2}).make();   //make DP distincive so that we can test the precendence properly
        final Genotype gGG = new GenotypeBuilder("4", GG).DP(10).AD(new int[]{1,9}).make();

        final double log10PError = -5;
        final double qual = -10.0 * log10PError;
        tests.add(new Object[]{new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make(), qual / 10});   //het
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
        final String QD = (String)annotatedMap.get(GATKVCFConstants.QUAL_BY_DEPTH_KEY);
        Assert.assertEquals(Double.valueOf(QD), expectedQD, 0.0001);

        Assert.assertEquals(new QualByDepth().getKeyNames(), Collections.singletonList(GATKVCFConstants.QUAL_BY_DEPTH_KEY));
        Assert.assertEquals(new QualByDepth().getDescriptions(), Collections.singletonList(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.QUAL_BY_DEPTH_KEY)));
    }

    @Test
    public void testUsingDP() {
        final Allele A = Allele.create("A", true);
        final Allele C = Allele.create("C");

        final List<Allele> AC = Arrays.asList(A, C);
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
        final PerReadAlleleLikelihoodMap map= new PerReadAlleleLikelihoodMap();

        final Allele A = Allele.create("A", true);
        final Allele C = Allele.create("C");
        final Allele G = Allele.create("G");

        final List<Allele> AC = Arrays.asList(A, C);
        final int readDepth = 20;
        final String sample1 = "sample1";
        final int dpDepth = 30; //Note: using a different value on purpose so that we can check that reads are preferred over DP
        final Genotype gAC = new GenotypeBuilder(sample1, AC).DP(dpDepth).make();

        final double log10PError = -5;
        final double qual = -10.0 * log10PError;

        final int n1A= readDepth;
        for (int i = 0; i < n1A; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"), "n1A_" + i);
            read.setMappingQuality(20);
            map.add(read, A, -1.0);
            map.add(read, C, -100.0);  //try to fool it - add another likelihood to same read
            map.add(read, G, -1000.0);  //and a third one
        }

        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();
        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = Collections.singletonMap(sample1, map);
        final Map<String, Object> annotatedMap = new QualByDepth().annotate(null, vc, perReadAlleleLikelihoodMap);
        Assert.assertNotNull(annotatedMap, vc.toString());
        final String QD = (String)annotatedMap.get(GATKVCFConstants.QUAL_BY_DEPTH_KEY);

        final double expectedQD = qual/readDepth;
        Assert.assertEquals(Double.valueOf(QD), expectedQD, 0.0001);

        //Now we test that when AD is present, it trumps everything
        final Genotype gAC_withAD = new GenotypeBuilder("1", AC).DP(dpDepth).AD(new int[]{5,5}).make();
        final VariantContext vc_withAD = new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).genotypes(Arrays.asList(gAC_withAD)).make();
        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap_withAD = Collections.singletonMap(sample1, map);
        final Map<String, Object> annotatedMap_withAD = new QualByDepth().annotate(null, vc_withAD, perReadAlleleLikelihoodMap_withAD);
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
    public void testNull(final VariantContext vc){
        final Map<String, Object> annotatedMap = new QualByDepth().annotate(null, vc, null);
        Assert.assertNull(annotatedMap);
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

}
