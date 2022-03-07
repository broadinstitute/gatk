package org.broadinstitute.hellbender.tools.walkers.genotyper;

import com.google.common.base.Strings;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import static org.testng.Assert.assertEquals;

/**
 * Test {@link org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAlleleCounts}
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class GenotypeAlleleCountsUnitTest {

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testFirstError() {
        GenotypeAlleleCounts.first(-1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAlleleIndexAtError() {
        final GenotypeAlleleCounts first = GenotypeAlleleCounts.first(2);
        first.alleleIndexAt(-1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAlleleRankForError() {
        final GenotypeAlleleCounts first = GenotypeAlleleCounts.first(2);
        first.alleleRankFor(-1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAlleleCountAtError() {
        final GenotypeAlleleCounts first = GenotypeAlleleCounts.first(2);
        first.alleleCountAt(-1);
    }

    @Test(dataProvider = "ploidyData")
    public void testFirst(final int ploidy) {
        final GenotypeAlleleCounts subject = GenotypeAlleleCounts.first(ploidy);
        Assert.assertNotNull(subject);
        assertEquals(subject.ploidy(), ploidy);
        assertEquals(subject.distinctAlleleCount(), 1);
        assertEquals(subject.log10CombinationCount(), 0.0);
        assertEquals(subject.alleleCountAt(0), ploidy);
        assertEquals(subject.alleleCountFor(0), ploidy);
        assertEquals(subject.alleleRankFor(0), 0);
        assertEquals(subject.alleleRankFor(1), -2);
        Assert.assertTrue(subject.containsAllele(0));
        Assert.assertFalse(subject.containsAllele(1));
        assertEquals(subject.alleleIndexAt(0), 0);
        assertEquals(subject.maximumAlleleIndex(), 0);
        assertEquals(subject.minimumAlleleIndex(), 0);
        assertEquals(subject.compareTo(subject), 0);
        assertEquals(subject, subject);
        assertEquals(subject.index(), 0);
        assertEquals(subject.asAlleleList(testAlleles), Collections.nCopies(ploidy, testAlleles.get(0)));
        for (int maximumAlleleIndex = 0; maximumAlleleIndex <= MAXIMUM_ALLELE_INDEX; maximumAlleleIndex++) {
            final int[] expected = new int[maximumAlleleIndex + 1];
            expected[0] = ploidy;
        }

        Assert.assertNotNull(subject.toString());
        assertEquals(subject.toUnphasedGenotypeString(), ploidy == 1 ? "0" : Strings.repeat("0/", ploidy - 1) + "0");
    }

    @Test(dataProvider = "ploidyDataWithZero", dependsOnMethods = "testFirst")
    public void testNext(final int ploidy) {
        if (ploidy == 0) {
            testNextZeroPloidy();
        } else if (ploidy == 1) {
            testNextOnePloidy();
        } else {
            testPloidyTwoOrMore(ploidy);
        }
    }

    @Test(dataProvider = "ploidyDataWithZero", dependsOnMethods = "testNext")
    public void testIncrease(final int ploidy) {
        if (ploidy == 0) {
            testNextZeroPloidyIncrease();
        } else if (ploidy == 1) {
            testNextOnePloidyIncrease();
        } else {
            testPloidyTwoOrMoreIncrease(ploidy);
        }
    }

    private void testNextZeroPloidy() {
        final GenotypeAlleleCounts first = GenotypeAlleleCounts.first(0);
        final GenotypeAlleleCounts next = first.next();
        assertEquals(first, next);
        assertEquals(first.compareTo(next), 0);
        assertEquals(next.compareTo(first), 0);
        assertEquals(next.distinctAlleleCount(), 0);
        assertEquals(next.ploidy(), 0);
        assertEquals(next.index(), 0);
        assertEquals(next.asAlleleList(testAlleles), Collections.EMPTY_LIST);

        first.increase();
        assertEquals(first, next);
    }

    private void testNextOnePloidy() {
        final GenotypeAlleleCounts ploidy2 = GenotypeAlleleCounts.first(2);
        GenotypeAlleleCounts current = GenotypeAlleleCounts.first(1);

        while (!current.containsAllele(MAXIMUM_ALLELE_INDEX + 1)) {
            final GenotypeAlleleCounts next = current.next();
            assertEquals(next.log10CombinationCount(), 0.0);
            assertEquals(next.minimumAlleleIndex(), next.maximumAlleleIndex());
            assertEquals(next.minimumAlleleIndex(), current.minimumAlleleIndex() + 1);
            assertEquals(next.alleleCountAt(0), 1);
            assertEquals(next.alleleIndexAt(0), next.minimumAlleleIndex());
            assertEquals(next.alleleRankFor(next.minimumAlleleIndex()), 0);
            assertEquals(next.alleleRankFor(next.minimumAlleleIndex() + 1), -2);
            assertEquals(next.alleleCountFor(next.minimumAlleleIndex()), 1);
            assertEquals(next.alleleCountFor(next.minimumAlleleIndex() + 1), 0);
            assertEquals(next.ploidy(), 1);

            Assert.assertTrue(next.compareTo(current) > 0);
            Assert.assertTrue(current.compareTo(next) < 0);
            assertEquals(next.compareTo(next), 0);
            assertEquals(next, next);
            Assert.assertNotEquals(next, null);
            Assert.assertNotEquals(next.toString(), next);
            Assert.assertNotEquals(ploidy2, next);
            Assert.assertNotEquals(next, ploidy2);
            Assert.assertNotEquals(current, next);
            Assert.assertNotEquals(next.hashCode(), current.hashCode());
            Assert.assertNotEquals(next, current);

            assertEquals(next.index(), current.index() + 1);
            assertEquals(next.ploidy(), current.ploidy());

            assertEquals(next.asAlleleList(testAlleles), Collections.singletonList(testAlleles.get(next.maximumAlleleIndex())));

            current = next;
        }
    }

    private void testPloidyTwoOrMore(final int ploidy) {
        if (ploidy < 2) {
            throw new IllegalArgumentException();
        }

        GenotypeAlleleCounts current = GenotypeAlleleCounts.first(ploidy);

        while (true) {
            final GenotypeAlleleCounts next = current.next();
            if (next.containsAllele(MAXIMUM_ALLELE_INDEX + 1)) {
                break;
            }

            // test log10CombinationCount
            if (ploidy == 2) {
                assertEquals(next.log10CombinationCount(), next.distinctAlleleCount() == 2 ? Math.log10(2) : 0.0);
            } else if (ploidy == 3) {
                assertEquals(next.log10CombinationCount(),
                        next.distinctAlleleCount() == 3 ? Math.log10(6) : (next.distinctAlleleCount() == 2 ? Math.log10(6) - Math.log10(2) : 0.0));
            } else {
                if (next.distinctAlleleCount() == 1) {
                    assertEquals(next.log10CombinationCount(), 0.0);
                } else if (next.distinctAlleleCount() == ploidy) {
                    assertEquals(next.log10CombinationCount(), MathUtils.logToLog10(CombinatoricsUtils.factorialLog(ploidy)));
                }
            }

            //test forEach
            final List<Integer> alleleCountsAsList = new ArrayList<>(next.distinctAlleleCount() * 2);
            final Set<Integer> absentAlleles = new HashSet<>();
            next.forEachAlleleIndexAndCount((alleleIndex, alleleCount) -> {
                alleleCountsAsList.add(alleleIndex);
                alleleCountsAsList.add(alleleCount);
            });
            next.forEachAbsentAlleleIndex(absentAlleles::add, MAXIMUM_ALLELE_INDEX + 1);

            assertEquals(absentAlleles.size(), MAXIMUM_ALLELE_INDEX + 1 - next.distinctAlleleCount());
            next.forEachAlleleIndexAndCount((index, count) -> Assert.assertFalse(absentAlleles.contains(index)));

            if (current.distinctAlleleCount() == 1) {
                assertEquals(next.maximumAlleleIndex(), current.maximumAlleleIndex() + 1);
                assertEquals(next.distinctAlleleCount(), 2);
                assertEquals(next.minimumAlleleIndex(), 0);
            } else {
                assertEquals(next.maximumAlleleIndex(), current.maximumAlleleIndex());
                assertEquals(next.minimumAlleleIndex(), current.alleleCountAt(0) > 1 ? 0
                        : current.alleleCountAt(0) == 1 ? current.minimumAlleleIndex() + 1 : current.minimumAlleleIndex());
            }

            // Checking on 0's new count and current.minAllele + 1 alleles.
            assertEquals(next.alleleCountFor(0), current.alleleCountFor(current.minimumAlleleIndex()) - 1);
            assertEquals(next.alleleCountFor(current.minimumAlleleIndex() + 1),
                    current.alleleCountFor(current.minimumAlleleIndex() + 1) + 1);

            // Checks current.minAllele count
            assertEquals(next.alleleCountFor(current.minimumAlleleIndex()),
                    current.minimumAlleleIndex() == 0 ? current.alleleCountAt(0) - 1 : 0);

            int totalCountSum = 0;
            final int[] expectedAlleleCountsByIndex = new int[Math.max(MAXIMUM_ALLELE_INDEX, next.maximumAlleleIndex()) + 1];
            for (int i = 0; i < next.distinctAlleleCount(); i++) {
                final int count = next.alleleCountAt(i);
                final int index = next.alleleIndexAt(i);
                expectedAlleleCountsByIndex[index] = count;
                // Check consistency of alleleCountAt(x) and alleleCountFor(alleleIndexAt(x))
                assertEquals(next.alleleCountFor(index), count);
                totalCountSum += count;
                // Check on counts of, in theory, unaffected allele counts.
                if (index > current.minimumAlleleIndex() + 1) {
                    assertEquals(next.alleleCountFor(index), current.alleleCountFor(index));
                }
            }
            assertEquals(totalCountSum, ploidy);

            Assert.assertTrue(next.compareTo(current) > 0);
            Assert.assertTrue(current.compareTo(next) < 0);
            assertEquals(next.compareTo(next), 0);
            assertEquals(next, next);
            Assert.assertNotEquals(current, next);
            Assert.assertNotEquals(next, current);
            assertEquals(next.index(), current.index() + 1);
            assertEquals(next.ploidy(), ploidy);

            //Check asAlleleList.
            final List<Allele> expectedList = new ArrayList<>(ploidy);
            for (int i = 0; i < next.distinctAlleleCount(); i++) {
                for (int j = 0; j < next.alleleCountAt(i); j++) {
                    expectedList.add(testAlleles.get(next.alleleIndexAt(i)));
                }
            }
            assertEquals(next.asAlleleList(testAlleles), expectedList);

            current = next;
        }
    }

    private void testNextZeroPloidyIncrease() {
        final GenotypeAlleleCounts first = GenotypeAlleleCounts.first(0);
        final GenotypeAlleleCounts next = first.copy();
        next.increase();
        assertEquals(first, next);
        assertEquals(first.compareTo(next), 0);
        assertEquals(next.compareTo(first), 0);
        assertEquals(next.distinctAlleleCount(), 0);
        assertEquals(next.ploidy(), 0);
        assertEquals(next.index(), 0);
    }

    private void testNextOnePloidyIncrease() {
        final GenotypeAlleleCounts next = GenotypeAlleleCounts.first(1);

        while (!next.containsAllele(MAXIMUM_ALLELE_INDEX + 1)) {
            final GenotypeAlleleCounts current = next.copy();
            next.increase();
            assertEquals(next.minimumAlleleIndex(), next.maximumAlleleIndex());
            assertEquals(next.minimumAlleleIndex(), current.minimumAlleleIndex() + 1);
            assertEquals(next.alleleCountAt(0), 1);
            assertEquals(next.alleleIndexAt(0), next.minimumAlleleIndex());
            assertEquals(next.alleleRankFor(next.minimumAlleleIndex()), 0);
            assertEquals(next.alleleRankFor(next.minimumAlleleIndex() + 1), -2);
            assertEquals(next.alleleCountFor(next.minimumAlleleIndex()), 1);
            assertEquals(next.alleleCountFor(next.minimumAlleleIndex() + 1), 0);
            assertEquals(next.ploidy(), 1);

            Assert.assertTrue(next.compareTo(current) > 0);
            Assert.assertTrue(current.compareTo(next) < 0);
            assertEquals(next.compareTo(next), 0);
            assertEquals(next, next);
            Assert.assertNotEquals(current, next);
            Assert.assertNotEquals(next, current);

            assertEquals(next.index(), current.index() + 1);
            assertEquals(next.ploidy(), current.ploidy());
        }
    }

    private void testPloidyTwoOrMoreIncrease(final int ploidy) {
        if (ploidy < 2) {
            throw new IllegalArgumentException();
        }

        final GenotypeAlleleCounts next = GenotypeAlleleCounts.first(ploidy);

        while (!next.containsAllele(MAXIMUM_ALLELE_INDEX + 1)) {
            final GenotypeAlleleCounts current = next.copy();
            next.increase();
            if (current.distinctAlleleCount() == 1) {
                assertEquals(next.maximumAlleleIndex(), current.maximumAlleleIndex() + 1);
                assertEquals(next.distinctAlleleCount(), 2);
                assertEquals(next.minimumAlleleIndex(), 0);
            } else {
                assertEquals(next.maximumAlleleIndex(), current.maximumAlleleIndex());
                assertEquals(next.minimumAlleleIndex(), current.alleleCountAt(0) > 1 ? 0
                        : current.alleleCountAt(0) == 1 ? current.minimumAlleleIndex() + 1 : current.minimumAlleleIndex());
            }

            // Checking on 0's new count and current.minAllele + 1 alleles.
            assertEquals(next.alleleCountFor(0), current.alleleCountFor(current.minimumAlleleIndex()) - 1);
            assertEquals(next.alleleCountFor(current.minimumAlleleIndex() + 1),
                    current.alleleCountFor(current.minimumAlleleIndex() + 1) + 1);

            // Checks current.minAllele count
            assertEquals(next.alleleCountFor(current.minimumAlleleIndex()),
                    current.minimumAlleleIndex() == 0 ? current.alleleCountAt(0) - 1 : 0);

            int totalCountSum = 0;
            final int[] expectedAlleleCountsByIndex = new int[Math.max(MAXIMUM_ALLELE_INDEX, next.maximumAlleleIndex()) + 1];
            for (int i = 0; i < next.distinctAlleleCount(); i++) {
                final int count = next.alleleCountAt(i);
                final int index = next.alleleIndexAt(i);
                expectedAlleleCountsByIndex[index] = count;
                // Check consistency of alleleCountAt(x) and alleleCountFor(alleleIndexAt(x))
                assertEquals(next.alleleCountFor(index), count);
                totalCountSum += count;
                // Check on counts of, in theory, unaffected allele counts.
                if (index > current.minimumAlleleIndex() + 1) {
                    assertEquals(next.alleleCountFor(index), current.alleleCountFor(index));
                }
            }
            assertEquals(totalCountSum, ploidy);

            Assert.assertTrue(next.compareTo(current) > 0);
            Assert.assertTrue(current.compareTo(next) < 0);
            assertEquals(next.compareTo(next), 0);
            assertEquals(next, next);
            Assert.assertNotEquals(current, next);
            Assert.assertNotEquals(next, current);
            assertEquals(next.index(), current.index() + 1);
            assertEquals(next.ploidy(), ploidy);
        }
    }

    private static final int MAXIMUM_ALLELE_INDEX = 10;

    private static final List<Allele> testAlleles;

    static {
        final StringBuilder sb = new StringBuilder(51);
        testAlleles = new ArrayList<>(51);
        sb.append('A');
        for (int i = 0; i <= 50; i++) {
            testAlleles.add(Allele.create(sb.toString().getBytes(), i == 0));
            sb.append('A');
        }
    }

    private static final int[] PLOIDY = {1, 2, 3, 7, 10};
    private static final int[] PLOIDY_WITH_ZERO = {0, 1, 2, 3, 7, 10};

    @DataProvider(name = "ploidyData")
    public Object[][] ploidyData() {
        final Object[][] result = new Object[PLOIDY.length][];
        for (int i = 0; i < PLOIDY.length; i++) {
            result[i] = new Object[]{PLOIDY[i]};
        }
        return result;
    }

    @DataProvider(name = "ploidyDataWithZero")
    public Object[][] ploidyDataWithZero() {
        final Object[][] result = new Object[PLOIDY_WITH_ZERO.length][];
        for (int i = 0; i < PLOIDY_WITH_ZERO.length; i++) {
            result[i] = new Object[]{PLOIDY_WITH_ZERO[i]};
        }
        return result;
    }

    private Object[][] genotypeAlleleCountDataWithToStringOutput() {
        return new Object[][]{
                new Object[]{1, new String[]{"0", "1", "2", "3", "4"}},
                new Object[]{2, new String[]{"0/0", "0/1", "1/1", "0/2", "1/2", "2/2", "0/3", "1/3", "2/3", "3/3"}},
                new Object[]{3, new String[]{"0/0/0", "0/0/1", "0/1/1", "1/1/1", "0/0/2", "0/1/2", "1/1/2", "0/2/2", "1/2/2", "2/2/2"}}};
    }

    @DataProvider
    public Iterator<Object[]> genotypeAlleleCountDataWithToStringOutputIterable() {
        final List<Object[]> tests = new ArrayList<>();
        for (final Object[] test : genotypeAlleleCountDataWithToStringOutput()) {
            final Integer ploidy = (Integer) test[0];
            GenotypeAlleleCounts current = GenotypeAlleleCounts.first(ploidy);
            for (final String expectedGenotypeString : (String[]) test[1]) {
                tests.add(new Object[]{current, expectedGenotypeString});
                current = current.next();
            }
        }
        return tests.iterator();
    }

    @Test(dataProvider = "genotypeAlleleCountDataWithToStringOutputIterable")
    public void testToString(final GenotypeAlleleCounts gac, final String expectedtoString) {
        assertEquals(gac.toString(), expectedtoString);
    }
}
