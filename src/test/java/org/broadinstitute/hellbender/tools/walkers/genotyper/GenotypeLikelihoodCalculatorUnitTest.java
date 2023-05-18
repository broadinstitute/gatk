package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoodsUnitTester;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Tests {@link GenotypesCache}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class GenotypeLikelihoodCalculatorUnitTest {
    
    @Test(dataProvider = "ploidyAndMaximumAlleleAndReadCountsData")
    public void testLikelihoodCalculation(final int ploidy, final int alleleCount, final int[] readCount) {
        final AlleleLikelihoods<GATKRead, Allele> readLikelihoods = ReadLikelihoodsUnitTester.readLikelihoods(alleleCount, readCount);
        final int sampleCount = readCount.length;
        for (int s = 0; s < sampleCount ; s++) {
            final LikelihoodMatrix<GATKRead, Allele> sampleLikelihoods = readLikelihoods.sampleMatrix(s);
            final GenotypeLikelihoods genotypeLikelihoods = GenotypeLikelihoodCalculator.log10GenotypeLikelihoods(ploidy, sampleLikelihoods);
            final double[] genotypeLikelihoodsDoubles = genotypeLikelihoods.getAsVector();
            for (final GenotypeAlleleCounts gac : GenotypeAlleleCounts.iterable(ploidy, alleleCount)) {
                Assert.assertNotNull(genotypeLikelihoods);
                final double[] readGenotypeLikelihoods = new double[sampleLikelihoods.evidenceCount()];
                for (int r = 0; r < sampleLikelihoods.evidenceCount(); r++) {
                    final double[] compoments = new double[gac.distinctAlleleCount()];
                    for (int ar = 0; ar < gac.distinctAlleleCount(); ar++) {
                        final int a = gac.alleleIndexAt(ar);
                        final int aCount = gac.alleleCountAt(ar);
                        final double readLk = sampleLikelihoods.get(a, r);
                        compoments[ar] = readLk + Math.log10(aCount);
                    }
                    readGenotypeLikelihoods[r] = MathUtils.approximateLog10SumLog10(compoments) - Math.log10(ploidy);
                }
                final double genotypeLikelihood = MathUtils.sum(readGenotypeLikelihoods);
                Assert.assertEquals(genotypeLikelihoodsDoubles[gac.index()], genotypeLikelihood, 0.0001 * Math.abs(genotypeLikelihood));
            }
        }
    }


    private static final int[] MAXIMUM_ALLELE = { 1, 2, 5, 6};

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
}
