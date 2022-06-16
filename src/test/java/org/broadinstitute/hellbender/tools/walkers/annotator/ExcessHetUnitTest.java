package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.IntStream;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

public final class ExcessHetUnitTest extends GATKBaseTest {
    private static final double DELTA_PRECISION = 1E-6;
    private static final double DELTA_PRECISION_PARSE = 1E-3;
    private static final Allele A_REF = Allele.create("A", true);
    private static final Allele T = Allele.create("T");
    private static final Allele C = Allele.create("C");
    private static final int[] HET_PL_ARRAY = {240, 0, 240};
    private static final int[] HOM_REF_PL_ARRAY = {0, 60, 600};

    @Override
    public String getToolTestDataDir() {
        return toolsTestDir + "walkers/annotator/";
    }


    private static Genotype makeG(final String sample, final Allele a1, final Allele a2, final int... pls) {
        return new GenotypeBuilder(sample, Arrays.asList(a1, a2)).PL(pls).make();
    }

    private static VariantContext makeVC(final String source, final List<Allele> alleles, final Genotype... genotypes) {
        final int start = 10;
        final int stop = start; // alleles.contains(ATC) ? start + 3 : start;
        return new VariantContextBuilder(source, "1", start, stop, alleles)
                .genotypes(Arrays.asList(genotypes))
                .unfiltered()
                .make();
    }

    @Test
    public void testPositiveZeroPhredScore() {
        // we generate a large number of hom-ref genotypes to drive the p-value to 1 and the Phred score to 0
        final int numHomref = 100;
        final VariantContext vc = makeVC("2", Arrays.asList(A_REF, T),
                IntStream.range(0, numHomref + 1).boxed().map(i -> makeG("s" + i, A_REF, T, 0, 1000, 1000)).toArray(Genotype[]::new));

        final double result = ExcessHet.calculateEH(vc, vc.getGenotypes()).getValue();
        final double expected = 0.; // calculated using the python implementation at https://github.com/broadinstitute/gatk/issues/7392#issue-959501711
        Assert.assertEquals(result, expected, DELTA_PRECISION, "Pass");

        final Map<String, Object> annots = new ExcessHet().annotate(null, vc, null);
        Assert.assertEquals(annots.keySet(), Collections.singleton(GATKVCFConstants.EXCESS_HET_KEY), "annots");
        Assert.assertEquals(annots.values().size(), 1, "size");
        Assert.assertEquals(Double.parseDouble((String) annots.values().iterator().next()), expected, DELTA_PRECISION_PARSE, "het");
        Assert.assertEquals((String) annots.values().iterator().next(), "0.0000", "zero"); // guards against -0.0000 being returned
    }

    @Test
    public void testExcessHetForMultiallelicVCCompoundHets() {
        // make sure that compound gets (with no ref) don't add to het count
        final VariantContext vc = makeVC("1", Arrays.asList(A_REF, T, C),
                makeG("s1", A_REF, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s2", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s3", T, C, 7099, 2530, 7099, 3056, 0, 14931),
                makeG("s4", A_REF, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s5", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s6", A_REF, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s7", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s8", A_REF, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s9", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s10", A_REF, T, 2530, 0, 7099, 366, 3056, 14931));

        final double result = ExcessHet.calculateEH(vc, vc.getGenotypes()).getValue();
        final double expected = 2.8389324060524004; // calculated using the python implementation at https://github.com/broadinstitute/gatk/issues/7392#issue-959501711
        Assert.assertEquals(result, expected, DELTA_PRECISION, "Pass");
    }

    @Test
    public void testExcessHetForMultiallelicVCCompoundHetsRefAltFlip() {
        // make sure that compound gets (with no ref) don't add to het count
        final VariantContext vc = makeVC("1", Arrays.asList(A_REF, T, C),
                makeG("s1", A_REF, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s2", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s3", T, C, 7099, 2530, 7099, 3056, 0, 14931),
                makeG("s4", T, A_REF, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s5", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s6", A_REF, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s7", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s8", A_REF, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s9", T, T, 7099, 2530, 0, 7099, 366, 3056),
                makeG("s10", A_REF, T, 2530, 0, 7099, 366, 3056, 14931));

        final double result = ExcessHet.calculateEH(vc, vc.getGenotypes()).getValue();
        final double expected = 2.8389324060524004; // calculated using the python implementation at https://github.com/broadinstitute/gatk/issues/7392#issue-959501711
        Assert.assertEquals(result, expected, DELTA_PRECISION, "Pass");
    }

    @Test
    public void testExcessHetForMultiallelicVCDifferentAlts() {
        // make sure that hets with different alternate alleles all get counted
        final VariantContext vc = makeVC("2", Arrays.asList(A_REF, T, C),
                makeG("s1", A_REF, C, 4878, 1623, 11297, 0, 7970, 8847),
                makeG("s2", A_REF, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s3", A_REF, T, 3382, 0, 6364, 1817, 5867, 12246),
                makeG("s4", A_REF, T, 2488, 0, 9110, 3131, 9374, 12505),
                makeG("s5", A_REF, C, 4530, 2006, 18875, 0, 6847, 23949),
                makeG("s6", A_REF, T, 5325, 0, 18692, 389, 16014, 24570),
                makeG("s7", A_REF, T, 2936, 0, 29743, 499, 21979, 38630),
                makeG("s8", A_REF, T, 6902, 0, 8976, 45, 5844, 9061),
                makeG("s9", A_REF, T, 5732, 0, 10876, 6394, 11408, 17802),
                makeG("s10", A_REF, T, 2780, 0, 25045, 824, 23330, 30939));

        final double result = ExcessHet.calculateEH(vc, vc.getGenotypes()).getValue();
        final double expected = 22.562985944843152; // calculated using the python implementation at https://github.com/broadinstitute/gatk/issues/7392#issue-959501711
        Assert.assertEquals(result, expected, DELTA_PRECISION, "Pass");

        // test the annotate method
        final Map<String, Object> annots = new ExcessHet().annotate(null, vc, null);
        Assert.assertEquals(annots.keySet(), Collections.singleton(GATKVCFConstants.EXCESS_HET_KEY), "annots");
        Assert.assertEquals(annots.values().size(), 1, "size");
        Assert.assertEquals(Double.parseDouble((String) annots.values().iterator().next()), expected, DELTA_PRECISION_PARSE, "het");
    }

    @Test
    public void testFounderIDsAndPedigreeFile() {
        // make sure that hets with different alternate alleles all get counted
        final VariantContext vc = makeVC("2", Arrays.asList(A_REF, T, C),
                makeG("s1", A_REF, C, 4878, 1623, 11297, 0, 7970, 8847),
                makeG("s2", A_REF, T, 2530, 0, 7099, 366, 3056, 14931),
                makeG("s3", A_REF, T, 3382, 0, 6364, 1817, 5867, 12246),
                makeG("s4", A_REF, T, 2488, 0, 9110, 3131, 9374, 12505),
                makeG("s5", A_REF, C, 4530, 2006, 18875, 0, 6847, 23949),
                makeG("s6", A_REF, T, 5325, 0, 18692, 389, 16014, 24570),
                makeG("s7", A_REF, T, 2936, 0, 29743, 499, 21979, 38630),
                makeG("s8", A_REF, T, 6902, 0, 8976, 45, 5844, 9061),
                makeG("s9", A_REF, T, 5732, 0, 10876, 6394, 11408, 17802),
                makeG("s10", A_REF, T, 2780, 0, 25045, 824, 23330, 30939));

        final Set<String> founderIDs = new HashSet<>(Arrays.asList("s1", "s2", "s3", "s4", "s5"));

        final double result = ExcessHet.calculateEH(vc, vc.getGenotypes(founderIDs)).getValue();
        final double expected = 8.96250562461638; // calculated using the python implementation at https://github.com/broadinstitute/gatk/issues/7392#issue-959501711
        Assert.assertEquals(result, expected, DELTA_PRECISION, "Pass");

        // test the annotate method with FounderIDs
        final Map<String, Object> annots = new ExcessHet(founderIDs).annotate(null, vc, null);
        Assert.assertEquals(annots.keySet(), Collections.singleton(GATKVCFConstants.EXCESS_HET_KEY), "annots");
        Assert.assertEquals(annots.values().size(), 1, "size");
        Assert.assertEquals(Double.parseDouble((String) annots.values().iterator().next()), expected, DELTA_PRECISION_PARSE, "het");

        // test the annotate method with a Pedigree File
        final Map<String, Object> annotsPedigree = new ExcessHet(getTestFileGATKPath("testPedigree.ped")).annotate(null, vc, null);
        Assert.assertEquals(annotsPedigree.keySet(), Collections.singleton(GATKVCFConstants.EXCESS_HET_KEY), "annots");
        Assert.assertEquals(annotsPedigree.values().size(), 1, "size");
        Assert.assertEquals(Double.parseDouble((String) annotsPedigree.values().iterator().next()), expected, DELTA_PRECISION_PARSE, "het");
    }
    
    @Test
    public void testSingletonVsCommonAllele() {
        final List<Genotype> allGTs = new ArrayList<>();
        final int numHomRefGTs = 10000;
        for (int i = 0; i < numHomRefGTs; i++) {
            allGTs.add(makeG("ref" + i, A_REF, A_REF, HOM_REF_PL_ARRAY));
        }

        allGTs.add(makeG("het0", A_REF, T, HET_PL_ARRAY));
        int numHetGTs = 1;

        final VariantContext singleton = makeVC("singleton", Arrays.asList(A_REF, T), allGTs.toArray(new Genotype[allGTs.size()]));
        final double singletonValue = ExcessHet.calculateEH(singleton, singleton.getGenotypes()).getValue();

        final int targetNumHetGTs = 20;
        for (int i = numHetGTs; i < targetNumHetGTs; i++) {
            allGTs.add(makeG("het" + i, A_REF, T, HET_PL_ARRAY));
        }

        final VariantContext common = makeVC("common", Arrays.asList(A_REF, T), allGTs.toArray(new Genotype[allGTs.size()]));
        final double EHcommon = ExcessHet.calculateEH(common, common.getGenotypes()).getValue();

        Assert.assertTrue(Math.abs(singletonValue) < Math.abs(EHcommon), String.format("singleton=%f common=%f", singletonValue, EHcommon));
    }

    @Test
    public void testLargeCohorts() {
        final List<Genotype> allGTs = new ArrayList<>();
        final int numHomRefGTs = 1000000;
        for (int i = 0; i < numHomRefGTs; i++) {
            allGTs.add(makeG("ref" + i, A_REF, A_REF, HOM_REF_PL_ARRAY));
        }

        allGTs.add(makeG("het0", A_REF, T, HET_PL_ARRAY));
        int numHetGTs = 1;

        final VariantContext singleton = makeVC("singleton", Arrays.asList(A_REF, T), allGTs.toArray(new Genotype[allGTs.size()]));
        final double singletonValue = ExcessHet.calculateEH(singleton, singleton.getGenotypes()).getValue();

        for (int i = numHetGTs; i < 100; i++) {
            allGTs.add(makeG("het" + i, A_REF, T, HET_PL_ARRAY));
            numHetGTs++;
        }

        final VariantContext hundredton = makeVC("hundredton", Arrays.asList(A_REF, T), allGTs.toArray(new Genotype[allGTs.size()]));
        final double hundredtonValue = ExcessHet.calculateEH(hundredton, hundredton.getGenotypes()).getValue();

        Assert.assertTrue(Math.abs(singletonValue) < Math.abs(hundredtonValue), String.format("singleton=%f hundredton=%f", singletonValue, hundredtonValue));

        for (int i = numHetGTs; i < numHomRefGTs; i++)
            allGTs.add(makeG("het" + i, A_REF, T, HET_PL_ARRAY));

        final VariantContext common = makeVC("common", Arrays.asList(A_REF, T), allGTs.toArray(new Genotype[allGTs.size()]));
        final double commonValue = ExcessHet.calculateEH(common, common.getGenotypes()).getValue();

        Assert.assertTrue(Math.abs(hundredtonValue) < Math.abs(commonValue), String.format("hundredton=%f common=%f", hundredtonValue, commonValue));
    }

    @Test
    public void testAllHetsForLargeCohorts() {
        final int numGTs = 1000000;

        final List<Genotype> singletonGTs = new ArrayList<>();
        for (int i = 0; i < numGTs; i++) {
            singletonGTs.add(makeG("ref" + i, A_REF, A_REF, HOM_REF_PL_ARRAY));
        }

        singletonGTs.add(makeG("het0", A_REF, T, HET_PL_ARRAY));

        final VariantContext singleton = makeVC("singleton", Arrays.asList(A_REF, T), singletonGTs.toArray(new Genotype[singletonGTs.size()]));
        final double singletonValue = ExcessHet.calculateEH(singleton, singleton.getGenotypes()).getValue();

        final List<Genotype> allHetGTs = new ArrayList<>();
        for (int i = 0; i < numGTs; i++) {
            allHetGTs.add(makeG("het" + i, A_REF, T, HET_PL_ARRAY));
        }

        final VariantContext allHet = makeVC("allHet", Arrays.asList(A_REF, T), allHetGTs.toArray(new Genotype[allHetGTs.size()]));
        final double hetsValue = ExcessHet.calculateEH(allHet, allHet.getGenotypes()).getValue();

        Assert.assertTrue(Math.abs(singletonValue) < Math.abs(hetsValue), String.format("singleton=%f allHets=%f", singletonValue, hetsValue));

        // Since all hets is such an extreme case and the sample size is large here, we know that the p-value should be 0
        Assert.assertEquals(hetsValue, ExcessHet.PHRED_SCALED_MIN_P_VALUE, DELTA_PRECISION, "P-value of 0 should be phred scaled to " + ExcessHet.PHRED_SCALED_MIN_P_VALUE);
    }

    @DataProvider(name = "smallSets")
    public Object[][] counts() {
        return new Object[][]{
                {1, 0, 0, 1.},
                {1, 1, 0, 1.},
                // The following test cases are taken from Table 1 of Wigginton et al.
                // at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1199378/.
                // Results were calculated using the python implementation at
                // https://github.com/broadinstitute/gatk/issues/7392#issue-959501711,
                // which explicitly calculates and sums the likelihoods rather than using the recurrence relation,
                // and checked against the values in the table. See further comments in the issue.
                {5, 87, 8, 0.9999999998661845},
                {17, 81, 2, 0.9304235455950692},
                {21, 79, 0, 0.3096036754839548},
                // The following test cases were also generated using the python implementation.
                {707, 359, 9, 0.},                      // python implementation returns 4.5660915932231366E-67, which we take as 0.
                {804, 599, 70, 0.},                     // python implementation returns 4.239662903004714E-24, which we take as 0.
                {684, 559, 629, 1.},                    // python implementation returns 1.000000000002042, which we take as 1.
                {192, 835, 763, 1.},                    // python implementation returns 1.0000000000013074, which we take as 1.
                {723, 277, 754, 0.9999982776347127},    // for testing to DELTA_PRECISION = 1E-6
                {537, 845, 72, 0.14763417247333427},
                {847, 431, 448, 0.7957628300988617},
                {659, 147, 910, 0.9664916975760656},
                {633, 938, 84, 0.04808183902489469},
                {881, 690, 292, 0.6731260236412528},
                {554, 371, 184, 0.1930042049715512},
                {552, 438, 207, 0.9367121245319978},
                {886, 812, 216, 0.1488013326311986},
                {537, 845, 72, 0.14763417247333427},
                {8700, 3186, 5918, 0.46111947799744585},
                // The following test cases were also generated using the python implementation and
                // further guard against the possibility of overflow as seen in
                // https://github.com/samtools/htsjdk/issues/44.
                // This implementation doesn't have a problem, since we cast to double when calculating the midpoint,
                // but see e.g. https://github.com/samtools/htsjdk/pull/70.
                {1998, 1, 998001, 0.7359430663641368},
                {19800, 100, 980100, 0.5267641864418983},
                {180000, 10000, 810000, 0.5028886002589862},
                {500000, 250000, 250000, 0.5009973527964665}
        };
    }


    @Test(dataProvider = "smallSets")
    public void smallSets(final int hetCount, final int homrefCount, final int homvarCount, final double expected) {
        final double actual = ExcessHet.exactTest(hetCount, homrefCount, homvarCount);
        Assert.assertEquals(actual, expected, DELTA_PRECISION, "Pass");
        // we also test swapping homvarCount <-> homrefCount, which should be symmetric
        final double actualSymmetric = ExcessHet.exactTest(hetCount, homvarCount, homrefCount);
        Assert.assertEquals(actualSymmetric, expected, DELTA_PRECISION, "Pass");
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
    public void testEmptyIfNoGenotypes() {
        final ExcessHet ann = new ExcessHet();
        final Map<String, Object> annotate = ann.annotate(null, when(mock(VariantContext.class).getGenotypesOrderedByName()).thenReturn(Collections.<Genotype>emptyList()).getMock(), null);
        Assert.assertTrue(annotate.isEmpty());
    }
}
