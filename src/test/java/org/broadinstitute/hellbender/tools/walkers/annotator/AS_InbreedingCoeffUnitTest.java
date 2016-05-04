package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

import java.util.*;

public final class AS_InbreedingCoeffUnitTest extends BaseTest{
    private static double DELTA_PRECISION = 0.001;
    private Allele Aref, T, C;
    private int[] hetPLs, homRefPLs;

    @BeforeSuite
    public void setup() {
        // alleles
        Aref = Allele.create("A", true);
        T = Allele.create("T");
        C = Allele.create("C");

        // simulating 20 reads with Q30 base qualities
        hetPLs = new int[] {240, 0, 240};
        homRefPLs = new int[] {0, 60, 600};
    }

    private Genotype makeG(String sample, Allele a1, Allele a2, int... pls) {
        return new GenotypeBuilder(sample, Arrays.asList(a1, a2)).PL(pls).make();
    }

    private VariantContext makeVC(String source, List<Allele> alleles, Genotype... genotypes) {
        int start = 10;
        int stop = start; // alleles.contains(ATC) ? start + 3 : start;
        return new VariantContextBuilder(source, "1", start, stop, alleles)
                .genotypes(Arrays.asList(genotypes))
                .filters((String)null)
                .make();
    }

    @Test
    public void testInbreedingCoeffForMultiallelicVC() {
        //make sure that compound gets (with no ref) don't add to het count
        VariantContext test1 = makeVC("1", Arrays.asList(Aref,T,C),
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

        final AS_InbreedingCoeff testClass1 = new AS_InbreedingCoeff();
        final double ICresult1 = testClass1.calculateIC(test1, T);
        Assert.assertEquals(ICresult1, -0.4285714, DELTA_PRECISION, "Pass");
        final double ICresult1b = testClass1.calculateIC(test1, C);
        Assert.assertEquals(ICresult1b, -0.05263, DELTA_PRECISION, "Pass");

        //make sure that hets with different alternate alleles all get counted
        VariantContext test2 = makeVC("2", Arrays.asList(Aref,T,C),
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

        final AS_InbreedingCoeff testClass2 = new AS_InbreedingCoeff();
        final double ICresult2 = testClass2.calculateIC(test2, T);
        Assert.assertEquals(ICresult2, -0.666666, DELTA_PRECISION, "Pass");
        final double ICresult2b = testClass2.calculateIC(test2, C);
        Assert.assertEquals(ICresult2b, -0.111129, DELTA_PRECISION, "Pass");

        //test the annotate method
        final Map<String, Object> annots = new AS_InbreedingCoeff().annotate(null, test2, null);
        Assert.assertEquals(annots.keySet(), Collections.singleton(GATKVCFConstants.AS_INBREEDING_COEFFICIENT_KEY), "annots");
        Assert.assertEquals(annots.values().size(), 1, "size");
        final String annotString = (String)annots.values().iterator().next();
        final String[] split = annotString.split(",");
        Assert.assertEquals(split.length, 2);
        Assert.assertEquals(Double.parseDouble(split[0]), -0.666666, DELTA_PRECISION, "ic");
        Assert.assertEquals(Double.parseDouble(split[1]), -0.111129, DELTA_PRECISION, "ic");

        if (false) { //TODO reenable when AS_InbreedingCoeff uses founders
            final Map<String, Object> annots3 = new AS_InbreedingCoeff(Collections.singleton("s1")).annotate(null, test2, null);
            Assert.assertNull(annots3);//not enough samples

            final Map<String, Object> annots4 = new AS_InbreedingCoeff(new HashSet<>(Arrays.asList("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9"))).annotate(null, test2, null);
            Assert.assertNull(annots4);//not enough samples
        }

    }

    @Test
    public void testSingletonVsCommonAllele() {

        final List<Genotype> allGTs = new ArrayList<>();
        final int numHomRefGTs = 10000;
        allGTs.add(makeG("het0", Aref, T, hetPLs));

        for ( int i = 0; i < numHomRefGTs; i++ )
            allGTs.add(makeG("ref" + i, Aref, Aref, homRefPLs));

        int numHetGTs = 1;

        final VariantContext singleton = makeVC("singleton", Arrays.asList(Aref, T), allGTs.toArray(new Genotype[allGTs.size()]));

        AS_InbreedingCoeff testClass = new AS_InbreedingCoeff();
        final double ICsingleton = testClass.calculateIC(singleton, T);

        final int targetNumHetGTs = 20;
        for ( int i = numHetGTs; i < targetNumHetGTs; i++ )
            allGTs.add(makeG("het" + i, Aref, T, hetPLs));

        final VariantContext common = makeVC("common", Arrays.asList(Aref, T), allGTs.toArray(new Genotype[allGTs.size()]));
        testClass = new AS_InbreedingCoeff();
        final double ICcommon = testClass.calculateIC(common, T);

        Assert.assertTrue(Math.abs(ICsingleton) < Math.abs(ICcommon), String.format("singleton=%f common=%f", ICsingleton, ICcommon));
    }

    @Test
    public void testLargeCohorts() {

        final List<Genotype> allGTs = new ArrayList<>();
        final int numHomRefGTs = 1000000;
        for ( int i = 0; i < numHomRefGTs; i++ )
            allGTs.add(makeG("ref" + i, Aref, Aref, homRefPLs));

        allGTs.add(makeG("het0", Aref, T, hetPLs));
        int numHetGTs = 1;

        final VariantContext singleton = makeVC("singleton", Arrays.asList(Aref, T), allGTs.toArray(new Genotype[allGTs.size()]));
        AS_InbreedingCoeff testClass = new AS_InbreedingCoeff();
        final double ICsingleton = testClass.calculateIC(singleton, T);

        for ( int i = numHetGTs; i < 100; i++ ) {
            allGTs.add(makeG("het" + i, Aref, T, hetPLs));
            numHetGTs++;
        }

        final VariantContext hundredton = makeVC("hundredton", Arrays.asList(Aref, T), allGTs.toArray(new Genotype[allGTs.size()]));
        testClass = new AS_InbreedingCoeff();
        final double IChundredton = testClass.calculateIC(hundredton, T);

        Assert.assertTrue(Math.abs(ICsingleton) < Math.abs(IChundredton), String.format("singleton=%f hundredton=%f", ICsingleton, IChundredton));

        for ( int i = numHetGTs; i < numHomRefGTs; i++ )
            allGTs.add(makeG("het" + i, Aref, T, hetPLs));

        final VariantContext common = makeVC("common", Arrays.asList(Aref, T), allGTs.toArray(new Genotype[allGTs.size()]));
        testClass = new AS_InbreedingCoeff();
        final double ICcommon = testClass.calculateIC(common, T);

        Assert.assertTrue(Math.abs(IChundredton) < Math.abs(ICcommon), String.format("hundredton=%f common=%f", IChundredton, ICcommon));
    }

    @Test
    public void testAllHetsForLargeCohorts() {

        final int numGTs = 1000000;

        final List<Genotype> singletonGTs = new ArrayList<>();
        for ( int i = 0; i < numGTs; i++ )
            singletonGTs.add(makeG("ref" + i, Aref, Aref, homRefPLs));

        singletonGTs.add(makeG("het0", Aref, T, hetPLs));

        final VariantContext singleton = makeVC("singleton", Arrays.asList(Aref, T), singletonGTs.toArray(new Genotype[singletonGTs.size()]));
        AS_InbreedingCoeff testClass = new AS_InbreedingCoeff();
        final double ICsingleton = testClass.calculateIC(singleton, T);

        final List<Genotype> allHetGTs = new ArrayList<>();
        for ( int i = 0; i < numGTs; i++ )
            allHetGTs.add(makeG("het" + i, Aref, T, hetPLs));

        final VariantContext allHet = makeVC("allHet", Arrays.asList(Aref, T), allHetGTs.toArray(new Genotype[allHetGTs.size()]));
        final double ICHets = testClass.calculateIC(allHet, T);

        Assert.assertTrue(Math.abs(ICsingleton) < Math.abs(ICHets), String.format("singleton=%f allHets=%f", ICsingleton, ICHets));
    }


    @Test
    public void testTinyCohorts() {

        final List<Genotype> allGTs = new ArrayList<>();
        final int numHomRefGTs = AS_InbreedingCoeff.MIN_SAMPLES - 2;
        for (int i = 0; i < numHomRefGTs; i++) {
            allGTs.add(makeG("ref" + i, Aref, Aref, homRefPLs));
        }

        allGTs.add(makeG("het0", Aref, T, hetPLs));
        final VariantContext singleton = makeVC("singleton", Arrays.asList(Aref, T), allGTs.toArray(new Genotype[allGTs.size()]));

        //test the annotate method
        final Map<String, Object> annots = new AS_InbreedingCoeff().annotate(null, singleton, null);
        Assert.assertNull(annots); //not enough samples
    }

    @Test
    public void testLabels(){
        Assert.assertEquals(new AS_InbreedingCoeff().getKeyNames(), Collections.singletonList(GATKVCFConstants.AS_INBREEDING_COEFFICIENT_KEY));
        Assert.assertEquals(new AS_InbreedingCoeff().getDescriptions(), Collections.singletonList(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_INBREEDING_COEFFICIENT_KEY)));
    }

}

