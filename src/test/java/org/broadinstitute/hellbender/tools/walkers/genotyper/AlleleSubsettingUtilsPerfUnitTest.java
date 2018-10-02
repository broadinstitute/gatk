package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class AlleleSubsettingUtilsPerfUnitTest extends GATKBaseTest {

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

        int numIterations = 1000000;
        long start = System.nanoTime();
        for (int i = 0; i < numIterations; i++) {

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
        long end = System.nanoTime();
        long totalTimeMillis = (end - start) / 1000000;
        System.out.println("Total time millis for " + numIterations + " iterations: " + totalTimeMillis);
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

        final Genotype base = new GenotypeBuilder("NA12878").DP(10).GQ(100).make();

        final double[] homG3AllelesPL = new double[]{-20, -10, -30, -40, -50, 0};  // AA, AC, CC, AG, CG, GG
        final int[] homG3AllelesAD = new int[]{0, 1, 21};  // AA, AC, CC, AG, CG, GG
        final int[] homG3AllelesSAC = new int[]{0, 0, 1, 1, 21, 21};  // AA, AC, CC, AG, CG, GG


        tests.add(new Object[]{
                new VariantContextBuilder(vcBase).alleles(ACG).genotypes(new GenotypeBuilder(base).alleles(AA).AD(homG3AllelesAD).PL(homG3AllelesPL).
                        attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, homG3AllelesSAC).make()).make(),
                new VariantContextBuilder(vcBase).alleles(AG).make(),
                Collections.singletonList(new GenotypeBuilder(base).alleles(AA).PL(new double[]{-20, -40, 0}).AD(new int[]{0, 21}).
                        attribute(GATKVCFConstants.STRAND_COUNT_BY_SAMPLE_KEY, new int[]{0, 0, 21, 21}).GQ(200).make())});

        return tests.toArray(new Object[][]{});
    }
}
