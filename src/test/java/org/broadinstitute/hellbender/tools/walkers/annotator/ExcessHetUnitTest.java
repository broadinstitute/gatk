package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

import static org.mockito.ArgumentMatchers.refEq;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

public final class ExcessHetUnitTest extends GATKBaseTest {
    private static double DELTA_PRECISION = .001;
    private Allele Aref= Allele.create("A", true);
    private Allele T = Allele.create("T");
    private Allele C = Allele.create("C");
    private int[] hetPLs = {240, 0, 240};
    private int[] homRefPLs= {0, 60, 600};

    @Override
    public String getToolTestDataDir() {
        return toolsTestDir + "walkers/annotator/";
    }


    private Genotype makeG(String sample, Allele a1, Allele a2, int... pls) {
        return new GenotypeBuilder(sample, Arrays.asList(a1, a2)).PL(pls).make();
    }

    private VariantContext makeVC(String source, List<Allele> alleles, Genotype... genotypes) {
        int start = 10;
        int stop = start; // alleles.contains(ATC) ? start + 3 : start;
        return new VariantContextBuilder(source, "1", start, stop, alleles)
                .genotypes(Arrays.asList(genotypes))
                .filters((String) null)
                .make();
    }

    @Test
    public void testExcessHetForMultiallelicVC_compondHets() {
        //make sure that compound gets (with no ref) don't add to het count
        VariantContext test1 = makeVC("1", Arrays.asList(Aref, T, C),
                makeG("s1", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s2", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s3", T, C, 7099, 2530, 7099, 3056, 0, 14931),
                makeG("s4", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s5", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s6", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s7", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s8", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s9", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s10", Aref, T, 2530, 0, 7099, 366, 3056, 14931));

        final double result = ExcessHet.calculateEH(test1, test1.getGenotypes()).getValue();
        Assert.assertEquals(result, 5.85, DELTA_PRECISION, "Pass");
    }

    @Test
    public void testExcessHetForMultiallelicVC_compondHetsRefAltFlip() {
        //make sure that compound gets (with no ref) don't add to het count
        VariantContext test1 = makeVC("1", Arrays.asList(Aref, T, C),
                makeG("s1", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s2", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s3", T, C, 7099, 2530, 7099, 3056, 0, 14931),
                makeG("s4", T, Aref, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s5", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s6", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s7", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s8", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s9", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s10", Aref, T, 2530, 0, 7099, 366, 3056, 14931));

        final double result = ExcessHet.calculateEH(test1, test1.getGenotypes()).getValue();
        Assert.assertEquals(result, 5.85, DELTA_PRECISION, "Pass");
    }

    @Test
    public void testExcessHetForMultiallelicVC_differentAlts() {
        //make sure that hets with different alternate alleles all get counted
        VariantContext test2 = makeVC("2", Arrays.asList(Aref, T, C),
                makeG("s1", Aref, C, 4878, 1623, 11297, 0, 7970, 8847),
                makeG("s2", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s3", Aref, T, 3382, 0, 6364, 1817, 5867, 12246),
                makeG("s4", Aref, T, 2488, 0, 9110, 3131, 9374, 12505),
                makeG("s5", Aref, C, 4530, 2006, 18875, 0, 6847, 23949),
                makeG("s6", Aref, T, 5325, 0, 18692, 389, 16014, 24570),
                makeG("s7", Aref, T, 2936, 0, 29743, 499, 21979, 38630),
                makeG("s8", Aref, T, 6902, 0, 8976, 45, 5844, 9061),
                makeG("s9", Aref, T, 5732, 0, 10876, 6394, 11408, 17802),
                makeG("s10", Aref, T, 2780, 0, 25045, 824, 23330, 30939));

        final double result2 = ExcessHet.calculateEH(test2, test2.getGenotypes()).getValue();
        final double result = 25.573;
        Assert.assertEquals(result2, result, DELTA_PRECISION, "Pass");

        //test the annotate method
        final Map<String, Object> annots = new ExcessHet().annotate(null, test2, null);
        Assert.assertEquals(annots.keySet(), Collections.singleton(GATKVCFConstants.EXCESS_HET_KEY), "annots");
        Assert.assertEquals(annots.values().size(), 1, "size");
        Assert.assertEquals(Double.parseDouble((String)annots.values().iterator().next()), result, DELTA_PRECISION, "het");
    }

    @Test
    public void testFounderIDsAndPedigreeFile() {
        //make sure that hets with different alternate alleles all get counted
        VariantContext test2 = makeVC("2", Arrays.asList(Aref, T, C),
                makeG("s1", Aref, C, 4878, 1623, 11297, 0, 7970, 8847),
                makeG("s2", Aref, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s3", Aref, T, 3382, 0, 6364, 1817, 5867, 12246),
                makeG("s4", Aref, T, 2488, 0, 9110, 3131, 9374, 12505),
                makeG("s5", Aref, C, 4530, 2006, 18875, 0, 6847, 23949),
                makeG("s6", Aref, T, 5325, 0, 18692, 389, 16014, 24570),
                makeG("s7", Aref, T, 2936, 0, 29743, 499, 21979, 38630),
                makeG("s8", Aref, T, 6902, 0, 8976, 45, 5844, 9061),
                makeG("s9", Aref, T, 5732, 0, 10876, 6394, 11408, 17802),
                makeG("s10", Aref, T, 2780, 0, 25045, 824, 23330, 30939));

        Set<String> founderIDs = new HashSet<String>();
        founderIDs.addAll(Arrays.asList("s1","s2","s3","s4","s5"));

        final double result2 = ExcessHet.calculateEH(test2, test2.getGenotypes(founderIDs)).getValue();
        final double result = 11.972;
        Assert.assertEquals(result2, result, DELTA_PRECISION, "Pass");

        //test the annotate method with FounderIDs
        Map<String, Object> annots = new ExcessHet(founderIDs).annotate(null, test2, null);
        Assert.assertEquals(annots.keySet(), Collections.singleton(GATKVCFConstants.EXCESS_HET_KEY), "annots");
        Assert.assertEquals(annots.values().size(), 1, "size");
        Assert.assertEquals(Double.parseDouble((String)annots.values().iterator().next()), result, DELTA_PRECISION, "het");

        //test the annotate method with a Pedigree File
        annots = new ExcessHet(getTestFile("testPedigree.ped")).annotate(null, test2, null);
        Assert.assertEquals(annots.keySet(), Collections.singleton(GATKVCFConstants.EXCESS_HET_KEY), "annots");
        Assert.assertEquals(annots.values().size(), 1, "size");
        Assert.assertEquals(Double.parseDouble((String)annots.values().iterator().next()), result, DELTA_PRECISION, "het");
    }
    
    @Test
    public void testSingletonVsCommonAllele() {

        final List<Genotype> allGTs = new ArrayList<>();
        final int numHomRefGTs = 10000;
        for (int i = 0; i < numHomRefGTs; i++) {
            allGTs.add(makeG("ref" + i, Aref, Aref, homRefPLs));
        }

        allGTs.add(makeG("het0", Aref, T, hetPLs));
        int numHetGTs = 1;

        final VariantContext singleton = makeVC("singleton", Arrays.asList(Aref, T), allGTs.toArray(new Genotype[allGTs.size()]));
        final double singletonValue = ExcessHet.calculateEH(singleton, singleton.getGenotypes()).getValue();

        final int targetNumHetGTs = 20;
        for (int i = numHetGTs; i < targetNumHetGTs; i++) {
            allGTs.add(makeG("het" + i, Aref, T, hetPLs));
        }

        final VariantContext common = makeVC("common", Arrays.asList(Aref, T), allGTs.toArray(new Genotype[allGTs.size()]));
        final double EHcommon = ExcessHet.calculateEH(common, common.getGenotypes()).getValue();

        Assert.assertTrue(Math.abs(singletonValue) < Math.abs(EHcommon), String.format("singleton=%f common=%f", singletonValue, EHcommon));
    }

    @Test
    public void testLargeCohorts() {

        final List<Genotype> allGTs = new ArrayList<>();
        final int numHomRefGTs = 1000000;
        for (int i = 0; i < numHomRefGTs; i++) {
            allGTs.add(makeG("ref" + i, Aref, Aref, homRefPLs));
        }

        allGTs.add(makeG("het0", Aref, T, hetPLs));
        int numHetGTs = 1;

        final VariantContext singleton = makeVC("singleton", Arrays.asList(Aref, T), allGTs.toArray(new Genotype[allGTs.size()]));
        final double singletonValue = ExcessHet.calculateEH(singleton, singleton.getGenotypes()).getValue();

        for (int i = numHetGTs; i < 100; i++) {
            allGTs.add(makeG("het" + i, Aref, T, hetPLs));
            numHetGTs++;
        }

        final VariantContext hundredton = makeVC("hundredton", Arrays.asList(Aref, T), allGTs.toArray(new Genotype[allGTs.size()]));
        final double hundredtonValue = ExcessHet.calculateEH(hundredton, hundredton.getGenotypes()).getValue();

        Assert.assertTrue(Math.abs(singletonValue) < Math.abs(hundredtonValue), String.format("singleton=%f hundredton=%f", singletonValue, hundredtonValue));

        for (int i = numHetGTs; i < numHomRefGTs; i++)
            allGTs.add(makeG("het" + i, Aref, T, hetPLs));

        final VariantContext common = makeVC("common", Arrays.asList(Aref, T), allGTs.toArray(new Genotype[allGTs.size()]));
        final double commonValue = ExcessHet.calculateEH(common, common.getGenotypes()).getValue();

        Assert.assertTrue(Math.abs(hundredtonValue) < Math.abs(commonValue), String.format("hundredton=%f common=%f", hundredtonValue, commonValue));
    }

    @Test
    public void testAllHetsForLargeCohorts() {

        final int numGTs = 1000000;

        final List<Genotype> singletonGTs = new ArrayList<>();
        for (int i = 0; i < numGTs; i++) {
            singletonGTs.add(makeG("ref" + i, Aref, Aref, homRefPLs));
        }

        singletonGTs.add(makeG("het0", Aref, T, hetPLs));

        final VariantContext singleton = makeVC("singleton", Arrays.asList(Aref, T), singletonGTs.toArray(new Genotype[singletonGTs.size()]));
        final double singletonValue = ExcessHet.calculateEH(singleton, singleton.getGenotypes()).getValue();

        final List<Genotype> allHetGTs = new ArrayList<>();
        for (int i = 0; i < numGTs; i++) {
            allHetGTs.add(makeG("het" + i, Aref, T, hetPLs));
        }

        final VariantContext allHet = makeVC("allHet", Arrays.asList(Aref, T), allHetGTs.toArray(new Genotype[allHetGTs.size()]));
        final double hetsValue = ExcessHet.calculateEH(allHet, allHet.getGenotypes()).getValue();

        Assert.assertTrue(Math.abs(singletonValue) < Math.abs(hetsValue), String.format("singleton=%f allHets=%f", singletonValue, hetsValue));

        //Since all hets is such an extreme case and the sample size is large here, we know that the p-value should be 0
        Assert.assertEquals(hetsValue, ExcessHet.PHRED_SCALED_MIN_P_VALUE, DELTA_PRECISION, String.format("P-value of 0 should be phred scaled to " + ExcessHet.PHRED_SCALED_MIN_P_VALUE));
    }

    @DataProvider(name = "smallSets")
    public Object[][] counts() {
        return new Object[][]{
                {1, 0, 0, 0.5},
                {1, 1, 0, 0.5},
                {1, 1, 1, 0.7},
                {4, 0, 0, 0.114},
                {2, 1, 1, 0.571},
                {0, 2, 2, 0.957},
                {1, 1, 40, 0.982},
                {3, 0, 39, 0.482},
        };
    }


    @Test(dataProvider = "smallSets")
    public void smallSets(int hetCount, int homrefCount, int homvarCount, double expected) {
        final double actual = ExcessHet.exactTest(hetCount, homrefCount, homvarCount);
        Assert.assertEquals(actual, expected, DELTA_PRECISION, "Pass");
    }

    @DataProvider(name = "illegalArgsForExactTest")
    public Object[][] illegalArgsForExactTest() {
        return new Object[][]{
                {-1, 1, 1},
                {1, -1, 1},
                {1, 1, -1},
        };
    }

    @Test(expectedExceptions = IllegalArgumentException.class, dataProvider = "illegalArgsForExactTest")
    public void testIllegalArgs(final int hetCount, final int refCount, final int homCount){
        ExcessHet.exactTest(hetCount, refCount, homCount);
    }

    @Test
    public void testLabels(){
        Assert.assertEquals(new ExcessHet().getKeyNames(), Collections.singletonList(GATKVCFConstants.EXCESS_HET_KEY));
        Assert.assertEquals(new ExcessHet().getDescriptions(), Collections.singletonList(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.EXCESS_HET_KEY)));
    }

    @Test
    public void testEmptyIfNoGenotypes() throws Exception {
        final ExcessHet ann = new ExcessHet();
        final Map<String, Object> annotate = ann.annotate(null, when(mock(VariantContext.class).getGenotypesOrderedByName()).thenReturn(Collections.<Genotype>emptyList()).getMock(), null);
        Assert.assertTrue(annotate.isEmpty());
    }

}
