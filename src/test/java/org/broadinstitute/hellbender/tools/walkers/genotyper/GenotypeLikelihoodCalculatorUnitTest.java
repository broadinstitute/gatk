package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoodsUnitTester;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Tests {@link GenotypeLikelihoodCalculators} and {@link GenotypeLikelihoodCalculator}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class GenotypeLikelihoodCalculatorUnitTest {

    @Test(dataProvider = "ploidyAndMaximumAlleleData")
    public void testPloidyAndMaximumAllele(final int ploidy, final int alleleCount) {
        final GenotypeLikelihoodCalculator calculator = new GenotypeLikelihoodCalculators().getInstance(ploidy, alleleCount);
        Assert.assertNotNull(calculator);
        Assert.assertEquals(calculator.ploidy(), ploidy);
        Assert.assertEquals(calculator.alleleCount(), alleleCount);
        Assert.assertEquals(calculator.genotypeCount(), calculateGenotypeCount(ploidy, alleleCount), " ploidy = " + ploidy + " alleleCount = " + alleleCount);
        final int genotypeCount = calculator.genotypeCount();
        final int testGenotypeCount = Math.min(30000, genotypeCount);
        for (int i = 0; i < testGenotypeCount; i++) {
            final GenotypeAlleleCounts alleleCounts = calculator.genotypeAlleleCountsAt(i);
            Assert.assertNotNull(alleleCounts);
            if (i > 0)
                Assert.assertTrue(calculator.genotypeAlleleCountsAt(i - 1).compareTo(alleleCounts) < 0);
            final int[] alleleArray = new int[ploidy];
            int index = 0;
            for (int j = 0; j < alleleCounts.distinctAlleleCount(); j++)
                Arrays.fill(alleleArray, index, index += alleleCounts.alleleCountAt(j), alleleCounts.alleleIndexAt(j));
            final int[] alleleCountArray = new int[alleleCounts.distinctAlleleCount() << 1];
            alleleCounts.copyAlleleCounts(alleleCountArray,0);
            Assert.assertEquals(index, ploidy);
            Assert.assertEquals(calculator.allelesToIndex(alleleArray), i);
            Assert.assertEquals(calculator.alleleCountsToIndex(alleleCountArray), i);
        }
    }

    @Test(dataProvider = "ploidyAndMaximumAlleleAndReadCountsData", dependsOnMethods = "testPloidyAndMaximumAllele")
    public void testLikelihoodCalculation(final int ploidy, final int alleleCount, final int[] readCount) {
        final ReadLikelihoods<Allele> readLikelihoods = ReadLikelihoodsUnitTester.readLikelihoods(alleleCount, readCount);
        final GenotypeLikelihoodCalculator calculator = new GenotypeLikelihoodCalculators().getInstance(ploidy, alleleCount);
        final int genotypeCount = calculator.genotypeCount();
        final int testGenotypeCount = Math.min(30000, genotypeCount);
        final int sampleCount = readCount.length;
        for (int s = 0; s < sampleCount ; s++) {
            final LikelihoodMatrix<Allele> sampleLikelihoods = readLikelihoods.sampleMatrix(s);
            final GenotypeLikelihoods genotypeLikelihoods = calculator.genotypeLikelihoods(sampleLikelihoods);
            final double[] genotypeLikelihoodsDoubles = genotypeLikelihoods.getAsVector();
            Assert.assertEquals(genotypeLikelihoodsDoubles.length, genotypeCount);
            for (int i = 0; i < testGenotypeCount; i++) {
                final GenotypeAlleleCounts genotypeAlleleCounts = calculator.genotypeAlleleCountsAt(i);
                Assert.assertNotNull(genotypeLikelihoods);
                final double[] readGenotypeLikelihoods = new double[sampleLikelihoods.numberOfReads()];
                for (int r = 0; r < sampleLikelihoods.numberOfReads(); r++) {
                    final double[] compoments = new double[genotypeAlleleCounts.distinctAlleleCount()];
                    for (int ar = 0; ar < genotypeAlleleCounts.distinctAlleleCount(); ar++) {
                        final int a = genotypeAlleleCounts.alleleIndexAt(ar);
                        final int aCount = genotypeAlleleCounts.alleleCountAt(ar);
                        final double readLk = sampleLikelihoods.get(a, r);
                        compoments[ar] = readLk + Math.log10(aCount);
                    }
                    readGenotypeLikelihoods[r] = MathUtils.approximateLog10SumLog10(compoments) - Math.log10(ploidy);
                }
                final double genotypeLikelihood = MathUtils.sum(readGenotypeLikelihoods);
                Assert.assertEquals(genotypeLikelihoodsDoubles[i], genotypeLikelihood, 0.0001);
            }
        }
    }

    @Test(dataProvider = "ploidyAndMaximumAlleleAndNewMaximumAlleleData")
    public void testGenotypeIndexMap(final int ploidy, final int oldAlleleCount, final int newAlleleCount) {
        final Random rnd = Utils.getRandomGenerator();
        final int maxAlleleCount = Math.max(oldAlleleCount, newAlleleCount);
        final int[] alleleMap = new int[newAlleleCount];
        final Map<Integer,Set<Integer>> reverseMap = new LinkedHashMap<>(oldAlleleCount);
        for (int i = 0; i < alleleMap.length; i++) {
            alleleMap[i] = rnd.nextInt(oldAlleleCount);
            if (reverseMap.get(alleleMap[i]) == null) reverseMap.put(alleleMap[i],new LinkedHashSet<>(6));
            reverseMap.get(alleleMap[i]).add(i);
        }
        final GenotypeLikelihoodCalculators calculators = new GenotypeLikelihoodCalculators();
        final GenotypeLikelihoodCalculator calculator = calculators.getInstance(ploidy, maxAlleleCount);

        final int[] genotypeIndexMap = calculator.genotypeIndexMap(alleleMap, calculators);
        Assert.assertNotNull(genotypeIndexMap);
        Assert.assertEquals(genotypeIndexMap.length, calculators.genotypeCount(ploidy, newAlleleCount));

        final GenotypeLikelihoodCalculator oldCalculator = calculators.getInstance(ploidy, oldAlleleCount);
        final GenotypeLikelihoodCalculator newCalculator = calculators.getInstance(ploidy, newAlleleCount);

        for (int i = 0; i < genotypeIndexMap.length; i++) {
            final GenotypeAlleleCounts oldCounts = oldCalculator.genotypeAlleleCountsAt(genotypeIndexMap[i]);
            final GenotypeAlleleCounts newCounts = newCalculator.genotypeAlleleCountsAt(i);
            final int[] reverseCounts = new int[oldAlleleCount];
            for (int j = 0; j < newCounts.distinctAlleleCount(); j++) {
                final int newIndex = newCounts.alleleIndexAt(j);
                final int newRepeats = newCounts.alleleCountAt(j);
                final int expectedOldIndex = alleleMap[newIndex];
                final int oldIndexRank = oldCounts.alleleRankFor(expectedOldIndex);
                Assert.assertNotEquals(oldIndexRank, -1);
                final int oldIndex = oldCounts.alleleIndexAt(oldIndexRank);
                final int oldRepeats = oldCounts.alleleCountAt(oldIndexRank);
                Assert.assertEquals(oldIndex, expectedOldIndex);
                // not necessarily the same count if two or more new alleles map the same old allele.
                Assert.assertTrue(oldRepeats >= newRepeats);
                reverseCounts[oldIndex] += newRepeats;
            }
            for (int j = 0; j < oldAlleleCount; j++)
                Assert.assertEquals(oldCounts.alleleCountFor(j), reverseCounts[j]);
        }
    }


    // Simple inefficient calculation of the genotype count given the ploidy.
    private int calculateGenotypeCount(final int ploidy, final int alleleCount) {
        if (ploidy == 0)
            return 0;
        else if (ploidy == 1)
            return alleleCount;
        else if (ploidy == 2)
            return ((alleleCount) * (alleleCount + 1)) >> 1;
        else if (alleleCount == 0)
            return 0;
        else {
            return calculateGenotypeCount(ploidy - 1, alleleCount) +
                        calculateGenotypeCount(ploidy, alleleCount - 1);
        }
    }

    private static final int[] MAXIMUM_ALLELE = { 1, 2, 5, 6 };

    private static final int[] PLOIDY = { 1, 2, 3, 20 };

    private static final int[][] READ_COUNTS = {
            { 10 , 100, 50 },
            { 0, 100, 10, 1 , 50 },
            { 1, 2, 3, 4, 20 },
            { 10, 0 },
    };

    @DataProvider(name="ploidyAndMaximumAlleleAndReadCountsData")
    public Object[][] ploidyAndMaximumAlleleAndReadCountsData() {
        final Object[][] result = new Object[PLOIDY.length * MAXIMUM_ALLELE.length * READ_COUNTS.length][];
        int index = 0;
        for (final int i : PLOIDY)
            for (final int j : MAXIMUM_ALLELE)
                for (final int[] k : READ_COUNTS)
                result[index++] = new Object[] { i, j, k };
        return result;
    }

    @DataProvider(name="ploidyAndMaximumAlleleData")
    public Object[][] ploidyAndMaximumAlleleData() {
        final Object[][] result = new Object[PLOIDY.length * MAXIMUM_ALLELE.length][];
        int index = 0;
        for (final int i : PLOIDY)
            for (final int j : MAXIMUM_ALLELE)
                    result[index++] = new Object[] { i, j };
        return result;
    }

    @DataProvider(name="ploidyAndMaximumAlleleAndNewMaximumAlleleData")
    public Object[][] ploidyAndMaximumAlleleAndNewMaximumAlleleData() {
        final List<Object[]> result = new ArrayList<>(PLOIDY.length * MAXIMUM_ALLELE.length * 20);
        for (final int i : PLOIDY)
            for (final int j : MAXIMUM_ALLELE)
                for (int k = 0; k < (i < 10? j * 2 : j + 1); k++)
                    result.add(new Object[] { i, j, k });
        return result.toArray(new Object[result.size()][]);
    }
}
