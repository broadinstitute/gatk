package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAlleleCounts;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

/**
 * Unit test for {@link MergedAlleleList}.
 */
public class MergedAlleleListUnitTest extends BaseTest {

    private AlleleList<Allele> allelesOf(final String str) {
        final String[] parts = str.split(",");
        final List<Allele> alleles = Arrays.stream(parts)
                .map(s -> Allele.create(s.replace("*", ""), s.contains("*")))
                .sorted(MergedAlleleList.ALLELE_COMPARATOR)
                .collect(Collectors.toList());
        return new IndexedAlleleList<>(alleles);
    }

    @Test(dataProvider = "mergedListBadData", expectedExceptions = IllegalArgumentException.class)
    public void testMergedBadDataList(final String leftStr, final String rightStr) {
        final AlleleList<Allele> left = allelesOf(leftStr);
        final AlleleList<Allele> right = allelesOf(rightStr);
        final MergedAlleleList<Allele> merged = MergedAlleleList.merge(left, right);
    }

    @Test(dataProvider = "mergedListData")
    public void testMergedList(final String leftStr, final String rightStr, final String exectedStr) {
        final AlleleList<Allele> left = allelesOf(leftStr);
        final AlleleList<Allele> right = allelesOf(rightStr);
        final AlleleList<Allele> expected = allelesOf(exectedStr);
        final MergedAlleleList<Allele> merged = MergedAlleleList.merge(left, right);
        Assert.assertEquals(merged.asListOfAlleles(), expected.asListOfAlleles());
    }

    @Test(dataProvider = "mergedListDataAndPloidy")
    public void testGenotypeNumberMapping(final String leftStr, final String rightStr, final String exectedStr, final int ploidy) {
        final AlleleList<Allele> left = allelesOf(leftStr);
        final AlleleList<Allele> right = allelesOf(rightStr);
        final AlleleList<Allele> expected = allelesOf(exectedStr);
        final MergedAlleleList<Allele> merged = MergedAlleleList.merge(left, right);
        final Allele mergedRefAllele = merged.getUniqueReferenceAllele();
        final Random rdn = new Random((((ploidy * 31) +left.hashCode() * 31) + right.hashCode()) * 31);
        final GenotypeLikelihoodCalculators calculators = new GenotypeLikelihoodCalculators();
        final GenotypeLikelihoodCalculator leftCalculator = calculators.getInstance(ploidy, left.numberOfAlleles());
        final GenotypeLikelihoodCalculator expectedCalculator = calculators.getInstance(ploidy, expected.numberOfAlleles());
        for (int i = 0; i < 10; i++) {
            final double[] likelihoods = rdn.doubles(leftCalculator.genotypeCount())
                    .map(Math::log10)
                    .toArray();
            normalizeLikelihoods(likelihoods);
            final double[] mappedLikelihoods = merged.mapGenotypeLikelihoods(left, ploidy, likelihoods);

            final boolean[] consumedLeftGenotype = new boolean[likelihoods.length];
            for (int g = 0; g < expectedCalculator.genotypeCount(); g++) {
                final GenotypeAlleleCounts mergedGACounts = expectedCalculator.genotypeAlleleCountsAt(g);
                final int[] leftAlleleIndexes = mergedGACounts.asAlleleList(merged.asListOfAlleles()).stream()
                        .mapToInt(a -> left.indexOfAllele(a, mergedRefAllele, true, true))
                        .filter(idx -> idx >= 0)
                        .toArray();
                if (leftAlleleIndexes.length == ploidy) {
                    final int leftGenotypeIndex = leftCalculator.allelesToIndex(leftAlleleIndexes);
                    consumedLeftGenotype[leftGenotypeIndex] = true;
                    Assert.assertEquals(mappedLikelihoods[g], likelihoods[leftGenotypeIndex]);
                }
            }
            for (int j = 0; j < consumedLeftGenotype.length; j++) {
                Assert.assertTrue(consumedLeftGenotype[j], "consumed j genotype");
            }
        }
    }

    @Test(dataProvider = "mergedListDataAndPloidy")
    public void testGenotypeNumberMappingRemovingAlleles(final String leftStr, final String rightStr, final String exectedStr, final int ploidy) {
        final AlleleList<Allele> left = allelesOf(leftStr);
        final AlleleList<Allele> right = allelesOf(rightStr);
        final AlleleList<Allele> expected = allelesOf(exectedStr);
        final Random rdn = new Random((((ploidy * 31) +left.hashCode() * 31) + right.hashCode()) * 31);
        final GenotypeLikelihoodCalculators calculators = new GenotypeLikelihoodCalculators();
        final GenotypeLikelihoodCalculator expectedCalculator = calculators.getInstance(ploidy, expected.numberOfAlleles());
        for (int i = 0; i < 3; i++) {
            final List<Allele> list = new ArrayList<>(right.asListOfAlleles());
            for (int j = 0; j < i; j++) {
                list.add(Allele.create("<SYMB_" + i + ">"));
            }
            final AlleleList<Allele> alternativeRight = new IndexedAlleleList<>(list);
            final MergedAlleleList<Allele> merged = MergedAlleleList.merge(left, alternativeRight);
            final Allele mergedRefAllele = merged.getUniqueReferenceAllele();
            final GenotypeLikelihoodCalculator leftCalculator = calculators.getInstance(ploidy, left.numberOfAlleles());
            for (int k = 0; k < 10; k++) {
                final double[] likelihoods = rdn.doubles(leftCalculator.genotypeCount())
                        .map(Math::log10)
                        .toArray();
                normalizeLikelihoods(likelihoods);
                final double[] mappedLikelihoods = merged.mapGenotypeLikelihoods(left, ploidy, likelihoods);

                final boolean[] consumedLeftGenotype = new boolean[likelihoods.length];
                for (int g = 0; g < expectedCalculator.genotypeCount(); g++) {
                    final GenotypeAlleleCounts mergedGACounts = expectedCalculator.genotypeAlleleCountsAt(g);
                    final int[] leftAlleleIndexes = mergedGACounts.asAlleleList(merged.asListOfAlleles()).stream()
                            .mapToInt(a -> left.indexOfAllele(a, mergedRefAllele, true, true))
                            .filter(idx -> idx >= 0)
                            .toArray();
                    if (leftAlleleIndexes.length == ploidy) {
                        final int leftGenotypeIndex = leftCalculator.allelesToIndex(leftAlleleIndexes);
                        consumedLeftGenotype[leftGenotypeIndex] = true;
                        Assert.assertEquals(mappedLikelihoods[g], likelihoods[leftGenotypeIndex]);
                    } else {
                        try {
                            Assert.assertTrue(mappedLikelihoods[g] < MathUtils.arrayMax(mappedLikelihoods));
                        } catch (final AssertionError ex) {
                            final double[] mappedLikelihoods2 = merged.mapGenotypeLikelihoods(left, ploidy, likelihoods);
                            throw ex;
                        }

                    }
                }
                for (int j = 0; j < consumedLeftGenotype.length; j++) {
                    final GenotypeAlleleCounts counts = expectedCalculator.genotypeAlleleCountsAt(j);
                    if (counts.asAlleleList(merged.asListOfAlleles()).stream()
                            .mapToInt(a -> left.indexOfAllele(a, mergedRefAllele, true, true))
                            .filter(idx -> idx >= 0)
                            .count() == ploidy) {
                        Assert.assertTrue(consumedLeftGenotype[j], "consumed j genotype");
                    }
                }
            }
        }

    }

    private void normalizeLikelihoods(final double[] likelihoods) {
        final double max = MathUtils.arrayMax(likelihoods);
        for (int i = 0; i < likelihoods.length; i++) {
            likelihoods[i] -= max;
        }
    }

    @Test(dataProvider = "mergedListData")
    public void testNumberOfAllelesAndGetAlleles(final String leftStr, final String rightStr, final String exectedStr) {
        final AlleleList<Allele> left = allelesOf(leftStr);
        final AlleleList<Allele> right = allelesOf(rightStr);
        final AlleleList<Allele> expected = allelesOf(exectedStr);
        final MergedAlleleList<Allele> merged = MergedAlleleList.merge(left, right);
        Assert.assertEquals(merged.numberOfAlleles(), expected.numberOfAlleles());
        for (int i = 0; i < merged.numberOfAlleles(); i++) {
            Assert.assertEquals(merged.getAllele(i), expected.getAllele(i));
        }
    }

    @Test(dataProvider = "mergedListData")
    public void testMapAlleleAnnotationIndexes(final String leftStr, final String rightStr, final String exectedStr) {
        final AlleleList<Allele> left = allelesOf(leftStr);
        final AlleleList<Allele> right = allelesOf(rightStr);
        final MergedAlleleList<Allele> merged = MergedAlleleList.merge(left, right);
        final Random rdn = new Random((leftStr.hashCode() *31 )+ rightStr.hashCode() * 31);
        for (int i = 0; i < 10; i++) {
            final int[] leftRandom = rdn.ints(left.numberOfAlleles()).toArray();
            final int[] leftMapped = merged.mapAlleleAnnotation(left, leftRandom, VCFHeaderLineCount.R,  -101);
            for (int j = 0; j < merged.numberOfAlleles(); j++) {
                final int leftAlleleIndex = left.indexOfAllele(merged.getAllele(j), true, true);
                final int leftValue = leftAlleleIndex < 0 ? -101 : leftRandom[leftAlleleIndex];
                final int resultValue = leftMapped[j];
                try {
                    Assert.assertEquals(leftValue, resultValue);
                } catch (final AssertionError err) {
                    final MergedAlleleList<Allele> merged2 = MergedAlleleList.merge(left, right);
                    final int[] leftMapped2 = merged.mapAlleleAnnotation(left, leftRandom, VCFHeaderLineCount.R, -101);
                    throw err;
                }
            }
        }
    }

    @Test(dataProvider = "mergedListData")
    public void testMapGenotypeAnnotationIndexes(final String leftStr, final String rightStr, final String exectedStr) {
        final AlleleList<Allele> left = allelesOf(leftStr);
        final AlleleList<Allele> right = allelesOf(rightStr);
        final MergedAlleleList<Allele> merged = MergedAlleleList.merge(left, right);
        final Random rdn = new Random((leftStr.hashCode() *31 )+ rightStr.hashCode() * 31);
        for (int i = 0; i < 10; i++) {
            final int[] leftRandom = rdn.ints(left.numberOfAlleles()).toArray();
            final int[] leftMapped = merged.mapAlleleAnnotation(left, leftRandom, VCFHeaderLineCount.R, -101);
            for (int j = 0; j < merged.numberOfAlleles(); j++) {
                final int leftAlleleIndex = left.indexOfAllele(merged.getAllele(j), true, true);
                final int leftValue = leftAlleleIndex < 0 ? -101 : leftRandom[leftAlleleIndex];
                final int resultValue = leftMapped[j];
                Assert.assertEquals(leftValue, resultValue);
            }
        }
    }

    @DataProvider
    public Object[][] mergedListData() {
        final List<Object[]> result = new ArrayList<>();

        result.add(new Object[] {"A*,C,G", "A*,T,TC", "A*,C,G,T,TC" });
        result.add(new Object[] {"ATTT*,ATT,A,CTT", "ATTTTT*,A", "ATTTTT*,A,ATT,ATTTT,CTTTT"});
        result.add(new Object[] {"ATTT*", "ATT*", "ATTT*"});
        result.add(new Object[] {"C*,A", "TA,CCT", "C*,A,CCT,TA"});
        result.add(new Object[] {"C", "A", "A,C"});
        result.add(new Object[] {"CC,T*,TA", "T*,A", "T*,A,CC,TA"});
        result.add(new Object[] {"CC,T*,TA", "T,A", "T*,A,CC,TA"});
        result.add(new Object[] {"C", "A,<NON_REF>", "A,C,<NON_REF>"});
        result.add(new Object[] {"C*,A", "TA,CCT,<NON_REF>", "C*,A,CCT,TA,<NON_REF>"});
        result.add(new Object[] {"C*,A,<NON_REF>", "TA,CCT,<NON_REF>", "C*,A,CCT,TA,<NON_REF>"});

        final List<Object[]> base = new ArrayList<>(result);

        for (final Object[] params : base) {
            result.add(new Object[] { params[1], params[0], params[2]});
        }

        for (final Object[] params : base) {
            result.add(new Object[] { params[0], params[0], params[0]});
        }

        for (final Object[] params : base) {
            result.add(new Object[] { params[0], params[0], params[0]});
        }

        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider
    public Object[][] mergedListDataAndPloidy() {
        final Object[][] mergedListData = mergedListData();
        final int[] ploidies = {1, 2, 3, 4, 5};
        final List<Object[]> result = new ArrayList<>(mergedListData.length * ploidies.length);
        for (final int ploidy : ploidies) {
            for (final Object[] mergedListDataCase : mergedListData) {
                final Object[] newCase = Arrays.copyOf(mergedListDataCase, mergedListDataCase.length + 1);
                newCase[mergedListDataCase.length] = ploidy;
                result.add(newCase);
            }
        }
        return result.toArray(new Object[result.size()][]);
    }

    @DataProvider
    public Object[][] mergedListBadData() {
        final List<Object[]> result = new ArrayList<>();

        result.add(new Object[] {"C*", "T*"});
        return result.toArray(new Object[result.size()][]);
    }
}
