package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class AFCalculationResultUnitTest extends GATKBaseTest {
    private static class MyTest {
        final double[] Ls, expectedPosteriors;

        private MyTest(double[] ls, double[] expectedPosteriors) {
            Ls = ls;
            this.expectedPosteriors = expectedPosteriors;
        }

        @Override
        public String toString() {
            return "Ls [" + Utils.join(",", Ls) + "] expectedPosteriors [" + Utils.join(",", expectedPosteriors) + "]";
        }
    }

    @DataProvider(name = "TestComputePosteriors")
    public Object[][] makeTestCombineGLs() {
        List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{new MyTest(log10Even, log10Even)});

        for ( double L0 = -1e9; L0 < 0.0; L0 /= 10.0 ) {
            for ( double L1 = -1e2; L1 < 0.0; L1 /= 100.0 ) {
                final double[] input = new double[]{L0, L1};
                final double[] expected = MathUtils.normalizeLog10(input);
                tests.add(new Object[]{new MyTest(input, expected)});
            }
        }

        for ( double bigBadL = -1e50; bigBadL < -1e200; bigBadL *= 10 ) {
            // test that a huge bad likelihood remains, even with a massive better result
            for ( final double betterL : Arrays.asList(-1000.0, -100.0, -10.0, -1.0, -0.1, -0.01, -0.001, 0.0)) {
                tests.add(new Object[]{new MyTest(new double[]{bigBadL, betterL}, new double[]{bigBadL, 0.0})});
                tests.add(new Object[]{new MyTest(new double[]{betterL, bigBadL}, new double[]{0.0, bigBadL})});
            }
        }

        // test that a modest bad likelihood with an ~0.0 value doesn't get lost
        for ( final double badL : Arrays.asList(-10000.0, -1000.0, -100.0, -10.0)) {
            tests.add(new Object[]{new MyTest(new double[]{badL, -1e-9}, new double[]{badL, 0.0})});
            tests.add(new Object[]{new MyTest(new double[]{-1e-9, badL}, new double[]{0.0, badL})});
        }

        // test that a non-ref site gets reasonable posteriors with an ~0.0 value doesn't get lost
        for ( final double nonRefL : Arrays.asList(-100.0, -50.0, -10.0, -9.0, -8.0, -7.0, -6.0, -5.0)) {
            tests.add(new Object[]{new MyTest(new double[]{0.0, nonRefL}, new double[]{0.0, nonRefL})});
        }

        return tests.toArray(new Object[][]{});
    }


    static final double[] log10Even = MathUtils.normalizeLog10(new double[]{0.5, 0.5});
    private static final Allele C = Allele.create("C");
    private static final Allele A = Allele.create("A", true);
    static final List<Allele> alleles = Arrays.asList(A, C);

    @Test(dataProvider = "TestComputePosteriors")
    private void testComputingPosteriors(final MyTest data) {
        final int[] zeroAC = {0};
        final AFCalculationResult result = new AFCalculationResult(zeroAC, alleles, data.Ls, log10Even, Collections.singletonMap(C, -1.0));

        Assert.assertEquals(result.getLog10PosteriorOfAFEq0(), data.expectedPosteriors[0], 1e-3, "AF = 0 not expected");
        Assert.assertEquals(result.getLog10PosteriorOfAFGT0(), data.expectedPosteriors[1], 1e-3, "AF > 0 not expected");

        Assert.assertEquals(result.getLog10PriorOfAFEq0(), log10Even[0], 1e-3, "prior for AF > 0 not expected");
        Assert.assertEquals(result.getLog10PriorOfAFGT0(), log10Even[1], 1e-3, "prior for AF > 0 not expected");

        Assert.assertEquals(result.getLog10LikelihoodOfAFEq0(), data.Ls[0], 1e-3, "likelihood for AF > 0 not expected");
        Assert.assertEquals(result.getLog10LikelihoodOfAFGT0(), data.Ls[1], 1e-3, "likelihood for AF > 0 not expected");

        Assert.assertEquals(result.getAllelesUsedInGenotyping(), alleles, "alleles are different");

        Assert.assertNotNull(result.toString());//just making sure it does not blow up, ignoring contents


        Assert.assertEquals(result.getAlleleCountAtMLE(C), zeroAC[0]);
        //getLog10PosteriorOfAFEq0ForAllele
        //withNewPriors

        Assert.assertEquals(result.getAlleleCountsOfMLE(), zeroAC, "getAlleleCountsOfMLE not as expected");
        final double[] actualPosteriors = {result.getLog10PosteriorOfAFEq0(), result.getLog10PosteriorOfAFGT0()};
        Assert.assertEquals(MathUtils.sumLog10(actualPosteriors), 1.0, 1e-3, "Posteriors don't sum to 1 with 1e-3 precision");
    }

    @DataProvider(name = "TestIsPolymorphic")
    public Object[][] makeTestIsPolymorphic() {
        List<Object[]> tests = new ArrayList<>();

        final List<Double> pValues = new LinkedList<>();
        for ( final double p : Arrays.asList(0.01, 0.1, 0.9, 0.99, 0.999, 1 - 1e-4, 1 - 1e-5, 1 - 1e-6) )
            for ( final double espilon : Arrays.asList(-1e-7, 0.0, 1e-7) )
                pValues.add(p + espilon);

        for ( final double pNonRef : pValues  ) {
            for ( final double pThreshold : pValues ) {
                final boolean shouldBePoly = pNonRef >= pThreshold;
                if ( pNonRef != pThreshold)
                    // let's not deal with numerical instability
                    tests.add(new Object[]{ pNonRef, pThreshold, shouldBePoly });
            }
        }

        return tests.toArray(new Object[][]{});
    }

    private AFCalculationResult makePolymorphicTestData(final double pNonRef) {
        return new AFCalculationResult(
                new int[]{0},
                alleles,
                MathUtils.normalizeLog10(new double[]{1 - pNonRef, pNonRef}),
                log10Even,
                Collections.singletonMap(C, Math.log10(1 - pNonRef)));
    }

    @Test(dataProvider = "TestIsPolymorphic")
    private void testIsPolymorphic(final double pNonRef, final double pThreshold, final boolean shouldBePoly) {
            final AFCalculationResult result = makePolymorphicTestData(pNonRef);
            final boolean actualIsPoly = result.isPolymorphic(C, Math.log10(1 - pThreshold));
            Assert.assertEquals(actualIsPoly, shouldBePoly,
                    "isPolymorphic with pNonRef " + pNonRef + " and threshold " + pThreshold + " returned "
                            + actualIsPoly + " but the expected result is " + shouldBePoly);
    }

    @Test(dataProvider = "TestIsPolymorphic")
    private void testIsPolymorphicQual(final double pNonRef, final double pThreshold, final boolean shouldBePoly) {
        final AFCalculationResult result = makePolymorphicTestData(pNonRef);
        final double qual = QualityUtils.phredScaleCorrectRate(pThreshold);
        final boolean actualIsPoly = result.isPolymorphicPhredScaledQual(C, qual);
        Assert.assertEquals(actualIsPoly, shouldBePoly,
                "isPolymorphic with pNonRef " + pNonRef + " and threshold " + pThreshold + " returned "
                        + actualIsPoly + " but the expected result is " + shouldBePoly);
    }

    @Test(dataProvider = "TestComputePosteriors")
    private void test(final MyTest data) {
        final AFCalculationResult result = new AFCalculationResult(new int[]{0}, alleles, data.Ls, log10Even, Collections.singletonMap(C, -1.0));

        Assert.assertEquals(result.getLog10PosteriorOfAFEq0(), data.expectedPosteriors[0], 1e-3, "AF = 0 not expected");
        Assert.assertEquals(result.getLog10PosteriorOfAFGT0(), data.expectedPosteriors[1], 1e-3, "AF > 0 not expected");

        final double[] actualPosteriors = {result.getLog10PosteriorOfAFEq0(), result.getLog10PosteriorOfAFGT0()};
        Assert.assertEquals(MathUtils.sumLog10(actualPosteriors), 1.0, 1e-3, "Posteriors don't sum to 1 with 1e-3 precision");
    }
}
