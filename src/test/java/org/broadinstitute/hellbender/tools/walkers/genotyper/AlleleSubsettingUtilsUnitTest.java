package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;


public class AlleleSubsettingUtilsUnitTest extends BaseTest {

    private static final Allele Aref = Allele.create("A", true);
    private static final Allele C = Allele.create("C");
    private static final Allele G = Allele.create("G");

    @Test(dataProvider = "updatePLsSACsAndADData")
    public void testUpdatePLsAndADData(final VariantContext originalVC,
                                       final VariantContext selectedVC,
                                       final List<Genotype> expectedGenotypes) {
        // initialize cache of allele anyploid indices
        for (final Genotype genotype : originalVC.getGenotypes()) {
            GenotypeLikelihoods.initializeAnyploidPLIndexToAlleleIndices(originalVC.getNAlleles() - 1, genotype.getPloidy());
        }

        final VariantContext selectedVCwithGTs = new VariantContextBuilder(selectedVC).genotypes(originalVC.getGenotypes()).make();

        final GenotypesContext oldGs = selectedVCwithGTs.getGenotypes();
        final GenotypesContext actual = selectedVCwithGTs.getNAlleles() == originalVC.getNAlleles() ? oldGs :
                                        AlleleSubsettingUtils.subsetAlleles(oldGs, 0, originalVC.getAlleles(),
                                                                            selectedVCwithGTs.getAlleles(),
                                                                            GenotypeAssignmentMethod.DO_NOT_ASSIGN_GENOTYPES,
                                                                            originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0));

        Assert.assertEquals(actual.size(), expectedGenotypes.size());
        for ( final Genotype expected : expectedGenotypes ) {
            final Genotype actualGT = actual.get(expected.getSampleName());
            Assert.assertNotNull(actualGT);
            VariantContextTestUtils.assertGenotypesAreEqual(actualGT, expected);
        }
    }

    @DataProvider(name = "updatePLsSACsAndADData")
    public Object[][] makeUpdatePLsSACsAndADData() {
        List<Object[]> tests = new ArrayList<>();

        final List<Allele> AA = Arrays.asList(Aref, Aref);
        final List<Allele> AC = Arrays.asList(Aref,C);
        final List<Allele> CC = Arrays.asList(C,C);
        final List<Allele> AG = Arrays.asList(Aref,G);
        final List<Allele> ACG = Arrays.asList(Aref,C,G);

        final VariantContext vcBase = new VariantContextBuilder("test", "20", 10, 10, AC).make();

        final double[] homRefPL = MathUtils.normalizeFromRealSpace(new double[]{0.9, 0.09, 0.01});
        final double[] hetPL = MathUtils.normalizeFromRealSpace(new double[]{0.09, 0.9, 0.01});
        final double[] homVarPL = MathUtils.normalizeFromRealSpace(new double[]{0.01, 0.09, 0.9});
        final double[] uninformative = new double[]{0, 0, 0};

        final Genotype base = new GenotypeBuilder("NA12878").DP(10).GQ(100).make();

        // the simple case where no selection occurs
        final Genotype aaGT = new GenotypeBuilder(base).alleles(AA).AD(new int[]{10,2}).PL(homRefPL).attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, new int[]{5, 10, 15, 20}).GQ(8).make();
        final Genotype acGT = new GenotypeBuilder(base).alleles(AC).AD(new int[]{10, 2}).PL(hetPL).attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, new int[]{5, 10, 15, 20}).GQ(8).make();
        final Genotype ccGT = new GenotypeBuilder(base).alleles(CC).AD(new int[]{10, 2}).PL(homVarPL).attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, new int[]{5, 10, 15, 20}).GQ(8).make();

        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(aaGT).make(), new VariantContextBuilder(vcBase).alleles(AC).make(), Collections.singletonList(new GenotypeBuilder(aaGT).make())});
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(acGT).make(), new VariantContextBuilder(vcBase).alleles(AC).make(), Collections.singletonList(new GenotypeBuilder(acGT).make())});
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(ccGT).make(), new VariantContextBuilder(vcBase).alleles(AC).make(), Collections.singletonList(new GenotypeBuilder(ccGT).make())});

        // uninformative test cases
        final Genotype uninformativeGT = new GenotypeBuilder(base).alleles(CC).noAD().PL(uninformative).GQ(0).make();
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(uninformativeGT).make(), new VariantContextBuilder(vcBase).alleles(AC).make(), Collections.singletonList(uninformativeGT)});
        final Genotype emptyGT = new GenotypeBuilder(base).alleles(GATKVariantContextUtils.noCallAlleles(2)).noAD().noPL().noGQ().make();
        tests.add(new Object[]{new VariantContextBuilder(vcBase).genotypes(emptyGT).make(), new VariantContextBuilder(vcBase).alleles(AC).make(), Collections.singletonList(emptyGT)});

        // actually subsetting down from multiple alt values
        final double[] homRef3AllelesPL = new double[]{0, -10, -20, -30, -40, -50};
        final double[] hetRefC3AllelesPL = new double[]{-10, 0, -20, -30, -40, -50};
        final double[] homC3AllelesPL = new double[]{-20, -10, 0, -30, -40, -50};
        final double[] hetRefG3AllelesPL = new double[]{-20, -10, -30, 0, -40, -50};
        final double[] hetCG3AllelesPL = new double[]{-20, -10, -30, -40, 0, -50}; // AA, AC, CC, AG, CG, GG
        final double[] homG3AllelesPL = new double[]{-20, -10, -30, -40, -50, 0};  // AA, AC, CC, AG, CG, GG

        final int[] homRef3AllelesAD = new int[]{20, 0, 1};
        final int[] hetRefC3AllelesAD = new int[]{10, 10, 1};
        final int[] homC3AllelesAD = new int[]{0, 20, 1};
        final int[] hetRefG3AllelesAD = new int[]{10, 0, 11};
        final int[] hetCG3AllelesAD = new int[]{0, 12, 11}; // AA, AC, CC, AG, CG, GG
        final int[] homG3AllelesAD = new int[]{0, 1, 21};  // AA, AC, CC, AG, CG, GG

        final double[] haploidRef3AllelesPL = new double[]{0, -10, -20};
        final double[] haploidAltC3AllelesPL = new double[]{-10, 0, -20};
        final double[] haploidAltG3AllelesPL = new double[]{-20, -10, 0};

        // for P=3 and N=2, the ordering is 000, 001, 011, 111, 002, 012, 112, 022, 122, 222
        final double[] triploidRef3AllelesPL = new double[]{0, -10, -20, -30, -40, -50, -60, -70, -80, -90};
        final double[] tripoidAltC3AllelesPL = new double[]{-10, 0, -20, -30, -40, -50, -60, -70, -80, -90};

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(Collections.singletonList(Aref)).AD(homRef3AllelesAD).PL(haploidRef3AllelesPL).make()).make(),
                new VariantContextBuilder(vcBase).alleles(AC).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(Collections.singletonList(Aref)).PL(new double[]{0, -10}).AD(new int[]{20, 0}).GQ(100).make())});
        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(Collections.singletonList(Aref)).AD(hetRefC3AllelesAD).PL(haploidAltC3AllelesPL).make()).make(),
                new VariantContextBuilder(vcBase).alleles(AC).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(Collections.singletonList(Aref)).PL(new double[]{-10, 0}).AD(new int[]{10, 10}).GQ(100).make())});
        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(Collections.singletonList(Aref)).AD(homC3AllelesAD).PL(haploidAltG3AllelesPL).make()).make(),
                new VariantContextBuilder(vcBase).alleles(AG).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(Collections.singletonList(Aref)).PL(new double[]{-20, 0}).AD(new int[]{0, 1}).GQ(200).make())});

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(Arrays.asList(Aref, Aref, Aref)).AD(homRef3AllelesAD).PL(triploidRef3AllelesPL).make()).make(),
                new VariantContextBuilder(vcBase).alleles(AC).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(Arrays.asList(Aref, Aref, Aref)).PL(new double[]{0, -10, -20, -30}).AD(new int[]{20, 0}).GQ(100).make())});

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(Arrays.asList(Aref, Aref, C)).AD(hetRefC3AllelesAD).PL(tripoidAltC3AllelesPL).make()).make(),
                new VariantContextBuilder(vcBase).alleles(AC).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(Arrays.asList(Aref, Aref, C)).PL(new double[]{-10, 0, -20, -30}).AD(new int[]{10, 10}).GQ(100).make())});

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(Arrays.asList(Aref, Aref, G)).AD(homRef3AllelesAD).PL(triploidRef3AllelesPL).make()).make(),
                new VariantContextBuilder(vcBase).alleles(AG).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(Arrays.asList(Aref, Aref, G)).PL(new double[]{0, -40, -70, -90}).AD(new int[]{20, 1}).GQ(400).make())});

        final int[] homRef3AllelesSAC = new int[]{20, 19, 0, 1, 3, 4};
        final int[] hetRefC3AllelesSAC = new int[]{10, 9, 10, 9, 1, 1};
        final int[] homC3AllelesSAC = new int[]{0, 0, 20, 20, 1, 1};
        final int[] hetRefG3AllelesSAC = new int[]{10, 10, 0, 0, 11, 11};
        final int[] hetCG3AllelesSAC = new int[]{0, 0, 12, 12, 11, 11}; // AA, AC, CC, AG, CG, GG
        final int[] homG3AllelesSAC = new int[]{0, 0, 1, 1, 21, 21};  // AA, AC, CC, AG, CG, GG

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).AD(homRef3AllelesAD).PL(homRef3AllelesPL).
                        attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, homRef3AllelesSAC).make()).make(),
                new VariantContextBuilder(vcBase).alleles(AC).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(AA).PL(new double[]{0, -10, -20}).AD(new int[]{20, 0}).
                        attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, new int[]{20, 19, 0, 1}).GQ(100).make())});

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).AD(hetRefC3AllelesAD).PL(hetRefC3AllelesPL).
                        attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, hetRefC3AllelesSAC).make()).make(),
                new VariantContextBuilder(vcBase).alleles(AC).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(AA).PL(new double[]{-10, 0, -20}).AD(new int[]{10, 10}).
                        attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, new int[]{10, 9, 10, 9}).GQ(100).make())});

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).AD(homC3AllelesAD).PL(homC3AllelesPL).
                        attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, homC3AllelesSAC).make()).make(),
                new VariantContextBuilder(vcBase).alleles(AC).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(AA).PL(new double[]{-20, -10, 0}).AD(new int[]{0, 20}).
                        attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, new int[]{0, 0, 20, 20}).GQ(100).make())});

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).AD(hetRefG3AllelesAD).PL(hetRefG3AllelesPL).
                        attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, hetRefG3AllelesSAC).make()).make(),
                new VariantContextBuilder(vcBase).alleles(AG).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(AA).PL(new double[]{-20, 0, -50}).AD(new int[]{10, 11}).
                        attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, new int[]{10, 10, 11, 11}).GQ(200).make())});

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).AD(hetCG3AllelesAD).PL(hetCG3AllelesPL).
                        attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, hetCG3AllelesSAC).make()).make(),
                new VariantContextBuilder(vcBase).alleles(AG).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(AA).PL(new double[]{0, -20, -30}).AD(new int[]{0, 11}).
                        attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, new int[]{0, 0, 11, 11}).GQ(200).make())});

        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).AD(homG3AllelesAD).PL(homG3AllelesPL).
                        attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, homG3AllelesSAC).make()).make(),
                new VariantContextBuilder(vcBase).alleles(AG).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(AA).PL(new double[]{-20, -40, 0}).AD(new int[]{0, 21}).
                        attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, new int[]{0, 0, 21, 21}).GQ(200).make())});

        return tests.toArray(new Object[][]{});
    }


    // regression test for https://github.com/broadinstitute/gatk/issues/2157
    @Test
    public void testCalculateMostLikelyAllelesTieDoesntRemoveAllTiedAlleles(){
        VariantContext vc = new VariantContextBuilder(null, "1", 100, 100, Arrays.asList(Aref, C, G))
                .genotypes(Arrays.asList(new GenotypeBuilder("sample1", Arrays.asList(C,G)).PL( new double[]{5, 5, 5, 5, 0, 5}).make())).make();
        Assert.assertEquals(AlleleSubsettingUtils.calculateMostLikelyAlleles(vc, 2, 1), Arrays.asList(Aref,C)) ;
    }

    @DataProvider
    public Object[][] getAllelesWithScores(){
        return new Object[][]{
                {1, Arrays.asList(Aref, C, G), new double[]{0,5,2}, Arrays.asList(Aref, C)}, //first is best
                {1, Arrays.asList(Aref, C, G), new double[]{0,2,5}, Arrays.asList(Aref, G)}, //second is best
                {1, Arrays.asList(Aref, C, G), new double[]{0,1,1}, Arrays.asList(Aref, C)}, //tie chooses first
                {1, Arrays.asList(Aref, C, G), new double[]{5,1,1}, Arrays.asList(Aref, C)}, //ref score is ignored chooses first
                {2, Arrays.asList(Aref, C, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE, G), new double[]{0,5,0,2}, Arrays.asList(Aref, C, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE, G) }, //keep NON_REF in order
                {1, Arrays.asList(Aref, C, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE, G), new double[]{0,5,0,2}, Arrays.asList(Aref, C, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE)}, //keep NON_REF in order when trimming
                {1, Arrays.asList(Aref, C, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE, G), new double[]{0,5,0,2}, Arrays.asList(Aref, C, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE)}, //keep NON_REF in order when trimming
        };
    }

    @Test(dataProvider = "getAllelesWithScores")
    public void testThatFilteringWorksCorrectly(int numToKeep, List<Allele> alleles, double[] scores, List<Allele> expected ){
        Assert.assertEquals(AlleleSubsettingUtils.filterToMaxNumberOfAltAllelesBasedOnScores(numToKeep, alleles, scores), expected);
    }

    @Test
    public void testCalculateMostLikelyAllelesPreconditions(){
        VariantContext vc = new VariantContextBuilder(null, "1", 100, 100, Arrays.asList(Aref, C, G)).make();
        Assert.assertThrows(IllegalArgumentException.class, () -> AlleleSubsettingUtils.calculateMostLikelyAlleles(null, 2, 2));
        Assert.assertThrows(IllegalArgumentException.class, () -> AlleleSubsettingUtils.calculateMostLikelyAlleles(vc, 0, 2));
        Assert.assertThrows(IllegalArgumentException.class, () -> AlleleSubsettingUtils.calculateMostLikelyAlleles(vc, 2, 0));
    }

    @Test
    public void testUninformativePLsAreKeptWhenDepthIsNotZero(){
        final Allele Aref = Allele.create("A", true);
        final List<Allele> alleles = Arrays.asList(Aref);
        final Genotype uniformativePL = new GenotypeBuilder("sample", alleles).PL(new int[] {0}).make();
        final GenotypesContext result  = AlleleSubsettingUtils.subsetAlleles(GenotypesContext.create(uniformativePL), 2,
                                                                      alleles, alleles, GenotypeAssignmentMethod.DO_NOT_ASSIGN_GENOTYPES, 10 );
        final Genotype genotype = result.get(0);
        Assert.assertTrue(genotype.hasPL());
        Assert.assertEquals(genotype.getPL(), new int[]{0});
    }

    @Test
    public void testCalculateLikelihoodSums() {
        // diploid, biallelic, three samples
        final List<Allele> twoAlleles = Arrays.asList(Aref, C);

        // the canonical ordering of genotypes is 0/0, 0/1, 1/1

        // sample 1 has GLs: {1.1, 0.1, 2.3}, so that its likeliest genotype has two copies of the alt allele
        // and its GL difference weight (see javadoc of the tested method) is 2.3 - 1.1 = 1.2

        // sample 2 will have GLs: {3.1, 0.1, 2.3}, so that its likeliest genotype is hom ref
        // and it contributes nothing to the likelihood sum

        // sample 3 will have GLs: {1.1, 4.1, 2.3}, so that its likeliest genotype has one copy of the alt allele
        // and its GL difference weight is 4.3 - 1.1 = 3.0

        // the total likelihood sum is thus 1.2 + 0.0 + 3.0 = 4.2

        final Genotype g1 = new GenotypeBuilder("sample1", twoAlleles).PL(new double[] {1.1, 0.1, 2.3}).make();
        final Genotype g2 = new GenotypeBuilder("sample2", twoAlleles).PL(new double[] {3.1, 0.1, 2.3}).make();
        final Genotype g3 = new GenotypeBuilder("sample3", twoAlleles).PL(new double[] {1.1, 4.1, 2.3}).make();
        final Genotype gNull = new GenotypeBuilder("sample4", twoAlleles).make();


        final VariantContext vc1 = new VariantContextBuilder("source", "contig", 1, 1, twoAlleles)
                .genotypes(Arrays.asList(g1, g2, g3, gNull)).make();

        Assert.assertEquals(AlleleSubsettingUtils.calculateLikelihoodSums(vc1, 2)[1], 4.2, 1.0e-8);

        // diploid, triallelic, two samples
        final List<Allele> threeAlleles = Arrays.asList(Aref, C, G);


        // the canonical ordering of genotypes is 0/0, 0/1, 1/1, 0/2, 1/2, 2/2

        // sample 1 has GLs: {0.0, 0.0, 0.0, 0.0, 3.1, 0.0}, so that its likeliest genotype has one copy of each
        // alt allele and its GL difference weight is 3.1
        // thus it contributes 3.1 to each alt allele

        // sample 2 has GLs: {0.0, 1.0, 0.0, 0.0, 0.0, 0.0}, so that its likeliest genotype has one copy of the first
        // alt allele and its GL difference weight is 1.0
        // thus it contributes 1.0 to the first alt allele

        // the likelihood sums are 4.1 and 3.1

        final Genotype g4 = new GenotypeBuilder("sample1", Arrays.asList(C, G)).PL(new double[] {0.0, 0.0, 0.0, 0.0, 3.1, 0.0}).make();
        final Genotype g5 = new GenotypeBuilder("sample2", Arrays.asList(Aref, C)).PL(new double[] {0.0, 1.0, 0.0, 0.0, 0.0, 0.0}).make();

        final VariantContext vc2 = new VariantContextBuilder("source", "contig", 1, 1, threeAlleles)
                .genotypes(Arrays.asList(g4, g5)).make();

        final double[] likelihoodSums2 = AlleleSubsettingUtils.calculateLikelihoodSums(vc2, 2);
        Assert.assertEquals(likelihoodSums2[1], 4.1, 1.0e-8);
        Assert.assertEquals(likelihoodSums2[2], 3.1, 1.0e-8);

        // triploid, biallelic, one sample
        // the canonical ordering of genotypes is 0/0/0, 0/0/1, 0/1/1, 1/1/1

        // sample 1 has GLs: {0.0, 0.0, 0.0, 3.5}, so that its likeliest genotype has three copies of the alt allele
        // and its GL difference weight is 3.5.

        final Genotype g6 = new GenotypeBuilder("sample1", Arrays.asList(C, C, C)).PL(new double[] {0.0, 0.0, 0.0, 3.5}).make();

        final VariantContext vc3 = new VariantContextBuilder("source", "contig", 1, 1, twoAlleles)
                .genotypes(Arrays.asList(g6)).make();

        Assert.assertEquals(AlleleSubsettingUtils.calculateLikelihoodSums(vc3, 3)[1], 3.5, 1.0e-8);
    }

}
