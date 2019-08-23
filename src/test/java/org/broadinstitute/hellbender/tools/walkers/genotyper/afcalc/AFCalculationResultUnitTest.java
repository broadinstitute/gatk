package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.util.*;

public final class AFCalculationResultUnitTest extends GATKBaseTest {


    private static final Allele A = Allele.create("A", true);
    private static final Allele C = Allele.create("C");
    private static final Allele T = Allele.create("T");
    private static final List<Allele> TWO_ALLELES = Arrays.asList(A, C);
    private static final List<Allele> THREE_ALLELES = Arrays.asList(A, C, T);


    @DataProvider(name = "AFCalculationResults")
    public Object[][] createAFCalculationResultTestData() {
        return new Object[][]{
                {new int[] {2}, TWO_ALLELES, 0.4, new double[] {0.4}},
                {new int[] {2,3}, THREE_ALLELES, 0.4, new double[] {0.5, 0.7}}
        };
    }

    @Test(dataProvider = "AFCalculationResults")
    private void test(final int[] mleCounts, final List<Allele> alleles, final double probabilityOfNoVariant, final double[] probabilityOfNoVariantByAllele) {
        final Map<Allele, Double> log10pRefByAllele = new HashMap<>();
        for (int n = 1; n < alleles.size(); n++) {
            log10pRefByAllele.put(alleles.get(n), Math.log10(probabilityOfNoVariantByAllele[n-1]));
        }
        final AFCalculationResult result = new AFCalculationResult(mleCounts, alleles, Math.log10(probabilityOfNoVariant), log10pRefByAllele);

        ArrayAsserts.assertArrayEquals(result.getAlleleCountsOfMLE(), mleCounts);
        Assert.assertEquals(result.getAllelesUsedInGenotyping(), alleles);

        for (int n = 1; n < alleles.size(); n++) {
            Assert.assertEquals(result.getAlleleCountAtMLE(alleles.get(n)), mleCounts[n-1]);
        }

        Assert.assertEquals(result.log10ProbVariantPresent(), Math.log10(1 - probabilityOfNoVariant), 1.0e-10);
    }
}