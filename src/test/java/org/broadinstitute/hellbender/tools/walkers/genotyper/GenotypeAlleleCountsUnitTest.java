package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Test {@link org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeAlleleCounts}
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class GenotypeAlleleCountsUnitTest {

    @Test(dataProvider="ploidyData")
    public void testFirst(final int ploidy) {
        final GenotypeAlleleCounts subject = GenotypeAlleleCounts.first(ploidy);
        Assert.assertNotNull(subject);
        Assert.assertEquals(subject.ploidy(), ploidy);
        Assert.assertEquals(subject.distinctAlleleCount(), 1);
        Assert.assertEquals(subject.alleleCountAt(0), ploidy);
        Assert.assertEquals(subject.alleleCountFor(0), ploidy);
        Assert.assertEquals(subject.alleleRankFor(0), 0);
        Assert.assertEquals(subject.alleleRankFor(1), -2);
        Assert.assertTrue(subject.containsAllele(0));
        Assert.assertFalse(subject.containsAllele(1));
        Assert.assertEquals(subject.alleleIndexAt(0), 0);
        Assert.assertEquals(subject.maximumAlleleIndex(), 0);
        Assert.assertEquals(subject.minimumAlleleIndex(), 0);
        Assert.assertTrue(subject.compareTo(subject) == 0);
        Assert.assertTrue(subject.equals(subject));
        Assert.assertEquals(subject.index(), 0);
        Assert.assertEquals(subject.asAlleleList(testAlleles), Collections.nCopies(ploidy, testAlleles.get(0)));
        for (int maximumAlleleIndex = 0; maximumAlleleIndex <= MAXIMUM_ALLELE_INDEX; maximumAlleleIndex++) {
            final int[] expected = new int[maximumAlleleIndex + 1];
            expected[0] = ploidy;
            Assert.assertEquals(subject.alleleCountsByIndex(maximumAlleleIndex), expected);
        }
    }

    @Test(dataProvider = "ploidyData",dependsOnMethods = "testFirst")
    public void testNext(final int ploidy) {
        if (ploidy == 0)
            testNextZeroPloidy();
        else if (ploidy == 1)
            testNextOnePloidy();
        else
            testPloidyTwoOrMore(ploidy);
    }

    @Test(dataProvider = "ploidyData",dependsOnMethods = "testNext")
    public void testIncrease(final int ploidy) {
        if (ploidy == 0)
            testNextZeroPloidyIncrease();
        else if (ploidy == 1)
            testNextOnePloidyIncrease();
        else
            testPloidyTwoOrMoreIncrease(ploidy);
    }

    private void testNextZeroPloidy() {
        final GenotypeAlleleCounts first = GenotypeAlleleCounts.first(0);
        final GenotypeAlleleCounts next = first.next();
        Assert.assertEquals(first, next);
        Assert.assertEquals(first.compareTo(next), 0);
        Assert.assertEquals(next.compareTo(first), 0);
        Assert.assertEquals(next.distinctAlleleCount(), 0);
        Assert.assertEquals(next.ploidy(), 0);
        Assert.assertEquals(next.index(), 0);
        Assert.assertEquals(next.asAlleleList(testAlleles), Collections.EMPTY_LIST);
        for (int maximumAlleleIndex = 0; maximumAlleleIndex <= 10; maximumAlleleIndex++) {
            final int[] expected = new int[maximumAlleleIndex + 1];
            Assert.assertEquals(next.alleleCountsByIndex(maximumAlleleIndex), expected);
        }
    }

    private void testNextOnePloidy() {
       final GenotypeAlleleCounts first = GenotypeAlleleCounts.first(1);
       GenotypeAlleleCounts current = first;

       while (!current.containsAllele(MAXIMUM_ALLELE_INDEX + 1)) {
            final GenotypeAlleleCounts next = current.next();
            Assert.assertEquals(next.minimumAlleleIndex(), next.maximumAlleleIndex());
            Assert.assertEquals(next.minimumAlleleIndex(), current.minimumAlleleIndex() + 1);
            Assert.assertEquals(next.alleleCountAt(0), 1);
            Assert.assertEquals(next.alleleIndexAt(0), next.minimumAlleleIndex());
            Assert.assertEquals(next.alleleRankFor(next.minimumAlleleIndex()), 0);
            Assert.assertEquals(next.alleleRankFor(next.minimumAlleleIndex() + 1), -2);
            Assert.assertEquals(next.alleleCountFor(next.minimumAlleleIndex()), 1);
            Assert.assertEquals(next.alleleCountFor(next.minimumAlleleIndex() + 1), 0);
            Assert.assertEquals(next.ploidy(), 1);

            Assert.assertTrue(next.compareTo(current) > 0);
            Assert.assertTrue(current.compareTo(next) < 0);
            Assert.assertTrue(next.compareTo(next) == 0);
            Assert.assertTrue(next.equals(next));
            Assert.assertFalse(next.equals(current));
            Assert.assertFalse(current.equals(next));

            Assert.assertEquals(next.index(), current.index() + 1);
            Assert.assertEquals(next.ploidy(), current.ploidy());

            Assert.assertEquals(next.asAlleleList(testAlleles), Collections.singletonList(testAlleles.get(next.maximumAlleleIndex())));

           for (int maximumAlleleIndex = 0; maximumAlleleIndex <= MAXIMUM_ALLELE_INDEX; maximumAlleleIndex++) {
               final int[] expected = new int[maximumAlleleIndex + 1];
               if (maximumAlleleIndex >= current.minimumAlleleIndex() + 1) expected[current.minimumAlleleIndex() + 1] = 1;
               Assert.assertEquals(next.alleleCountsByIndex(maximumAlleleIndex), expected);
           }
           current = next;
       }
    }

    private void testPloidyTwoOrMore(final int ploidy) {
        if (ploidy < 2)
            throw new IllegalArgumentException();

        GenotypeAlleleCounts current = GenotypeAlleleCounts.first(ploidy);

        while (!current.containsAllele(MAXIMUM_ALLELE_INDEX + 1)) {
            final GenotypeAlleleCounts next = current.next();
            if (current.distinctAlleleCount() == 1) {
                Assert.assertEquals(next.maximumAlleleIndex(), current.maximumAlleleIndex() + 1);
                Assert.assertEquals(next.distinctAlleleCount(), 2);
                Assert.assertEquals(next.minimumAlleleIndex(), 0);
            } else {
                Assert.assertEquals(next.maximumAlleleIndex(), current.maximumAlleleIndex());
                Assert.assertEquals(next.minimumAlleleIndex(), current.alleleCountAt(0) > 1 ? 0
                        : current.alleleCountAt(0) == 1 ? current.minimumAlleleIndex() + 1 : current.minimumAlleleIndex());
            }

            // Checking on 0's new count and current.minAllele + 1 alleles.
            Assert.assertEquals(next.alleleCountFor(0), current.alleleCountFor(current.minimumAlleleIndex()) - 1);
            Assert.assertEquals(next.alleleCountFor(current.minimumAlleleIndex() + 1),
                    current.alleleCountFor(current.minimumAlleleIndex() + 1) + 1);

            // Checks current.minAllele count
            Assert.assertEquals(next.alleleCountFor(current.minimumAlleleIndex()),
                    current.minimumAlleleIndex() == 0 ? current.alleleCountAt(0) - 1 : 0);

            int totalCountSum = 0;
            final int[] expectedAlleleCountsByIndex = new int[Math.max(MAXIMUM_ALLELE_INDEX, next.maximumAlleleIndex()) + 1];
            for (int i = 0; i < next.distinctAlleleCount(); i++) {
                final int count = next.alleleCountAt(i);
                final int index = next.alleleIndexAt(i);
                expectedAlleleCountsByIndex[index] = count;
                // Check consistency of alleleCountAt(x) and alleleCountFor(alleleIndexAt(x))
                Assert.assertEquals(next.alleleCountFor(index), count);
                totalCountSum += count;
                // Check on counts of, in theory, unaffected allele counts.
                if (index > current.minimumAlleleIndex() + 1)
                    Assert.assertEquals(next.alleleCountFor(index), current.alleleCountFor(index));
            }
            Assert.assertTrue(Arrays.equals(next.alleleCountsByIndex(Math.max(MAXIMUM_ALLELE_INDEX, next.maximumAlleleIndex())), expectedAlleleCountsByIndex));
            Assert.assertEquals(totalCountSum, ploidy);

            Assert.assertTrue(next.compareTo(current) > 0);
            Assert.assertTrue(current.compareTo(next) < 0);
            Assert.assertTrue(next.compareTo(next) == 0);
            Assert.assertTrue(next.equals(next));
            Assert.assertFalse(next.equals(current));
            Assert.assertFalse(current.equals(next));
            Assert.assertEquals(next.index(), current.index() + 1);
            Assert.assertEquals(next.ploidy(), ploidy);

            //Check asAlleleList.
            final List<Allele> expectedList = new ArrayList<>(ploidy);
            for (int i = 0; i < next.distinctAlleleCount(); i++) {
                for (int j = 0; j < next.alleleCountAt(i); j++) {
                    expectedList.add(testAlleles.get(next.alleleIndexAt(i)));
                }
            }
            Assert.assertEquals(next.asAlleleList(testAlleles), expectedList);

            current = next;
        }
    }

    private void testNextZeroPloidyIncrease() {
        final GenotypeAlleleCounts first = GenotypeAlleleCounts.first(0);
        final GenotypeAlleleCounts next = first.copy();
        next.increase();
        Assert.assertEquals(first, next);
        Assert.assertEquals(first.compareTo(next), 0);
        Assert.assertEquals(next.compareTo(first), 0);
        Assert.assertEquals(next.distinctAlleleCount(), 0);
        Assert.assertEquals(next.ploidy(), 0);
        Assert.assertEquals(next.index(), 0);
        for (int maximumAlleleIndex = 0; maximumAlleleIndex <= 10; maximumAlleleIndex++) {
            final int[] expected = new int[maximumAlleleIndex + 1];
            Assert.assertEquals(next.alleleCountsByIndex(maximumAlleleIndex), expected);
        }
    }

    private void testNextOnePloidyIncrease() {
        final GenotypeAlleleCounts first = GenotypeAlleleCounts.first(1);
        GenotypeAlleleCounts next = first;

        while (!next.containsAllele(MAXIMUM_ALLELE_INDEX + 1)) {
            final GenotypeAlleleCounts current = next.copy();
            next.increase();
            Assert.assertEquals(next.minimumAlleleIndex(), next.maximumAlleleIndex());
            Assert.assertEquals(next.minimumAlleleIndex(), current.minimumAlleleIndex() + 1);
            Assert.assertEquals(next.alleleCountAt(0), 1);
            Assert.assertEquals(next.alleleIndexAt(0), next.minimumAlleleIndex());
            Assert.assertEquals(next.alleleRankFor(next.minimumAlleleIndex()), 0);
            Assert.assertEquals(next.alleleRankFor(next.minimumAlleleIndex() + 1), -2);
            Assert.assertEquals(next.alleleCountFor(next.minimumAlleleIndex()), 1);
            Assert.assertEquals(next.alleleCountFor(next.minimumAlleleIndex() + 1), 0);
            Assert.assertEquals(next.ploidy(), 1);

            Assert.assertTrue(next.compareTo(current) > 0);
            Assert.assertTrue(current.compareTo(next) < 0);
            Assert.assertTrue(next.compareTo(next) == 0);
            Assert.assertTrue(next.equals(next));
            Assert.assertFalse(next.equals(current));
            Assert.assertFalse(current.equals(next));

            Assert.assertEquals(next.index(), current.index() + 1);
            Assert.assertEquals(next.ploidy(), current.ploidy());

            for (int maximumAlleleIndex = 0; maximumAlleleIndex <= MAXIMUM_ALLELE_INDEX; maximumAlleleIndex++) {
                final int[] expected = new int[maximumAlleleIndex + 1];
                if (maximumAlleleIndex >= current.minimumAlleleIndex() + 1) expected[current.minimumAlleleIndex() + 1] = 1;
                Assert.assertEquals(next.alleleCountsByIndex(maximumAlleleIndex), expected);
            }
        }
    }

    private void testPloidyTwoOrMoreIncrease(final int ploidy) {
        if (ploidy < 2)
            throw new IllegalArgumentException();

        GenotypeAlleleCounts next = GenotypeAlleleCounts.first(ploidy);

        while (!next.containsAllele(MAXIMUM_ALLELE_INDEX + 1)) {
            final GenotypeAlleleCounts current = next.copy();
            next.increase();
            if (current.distinctAlleleCount() == 1) {
                Assert.assertEquals(next.maximumAlleleIndex(), current.maximumAlleleIndex() + 1);
                Assert.assertEquals(next.distinctAlleleCount(), 2);
                Assert.assertEquals(next.minimumAlleleIndex(), 0);
            } else {
                Assert.assertEquals(next.maximumAlleleIndex(), current.maximumAlleleIndex());
                Assert.assertEquals(next.minimumAlleleIndex(), current.alleleCountAt(0) > 1 ? 0
                        : current.alleleCountAt(0) == 1 ? current.minimumAlleleIndex() + 1 : current.minimumAlleleIndex());
            }

            // Checking on 0's new count and current.minAllele + 1 alleles.
            Assert.assertEquals(next.alleleCountFor(0), current.alleleCountFor(current.minimumAlleleIndex()) - 1);
            Assert.assertEquals(next.alleleCountFor(current.minimumAlleleIndex() + 1),
                    current.alleleCountFor(current.minimumAlleleIndex() + 1) + 1);

            // Checks current.minAllele count
            Assert.assertEquals(next.alleleCountFor(current.minimumAlleleIndex()),
                    current.minimumAlleleIndex() == 0 ? current.alleleCountAt(0) - 1 : 0);

            int totalCountSum = 0;
            final int[] expectedAlleleCountsByIndex = new int[Math.max(MAXIMUM_ALLELE_INDEX, next.maximumAlleleIndex()) + 1];
            for (int i = 0; i < next.distinctAlleleCount(); i++) {
                final int count = next.alleleCountAt(i);
                final int index = next.alleleIndexAt(i);
                expectedAlleleCountsByIndex[index] = count;
                // Check consistency of alleleCountAt(x) and alleleCountFor(alleleIndexAt(x))
                Assert.assertEquals(next.alleleCountFor(index), count);
                totalCountSum += count;
                // Check on counts of, in theory, unaffected allele counts.
                if (index > current.minimumAlleleIndex() + 1)
                    Assert.assertEquals(next.alleleCountFor(index), current.alleleCountFor(index));
            }
            Assert.assertTrue(Arrays.equals(next.alleleCountsByIndex(Math.max(MAXIMUM_ALLELE_INDEX, next.maximumAlleleIndex())), expectedAlleleCountsByIndex));
            Assert.assertEquals(totalCountSum, ploidy);

            Assert.assertTrue(next.compareTo(current) > 0);
            Assert.assertTrue(current.compareTo(next) < 0);
            Assert.assertTrue(next.compareTo(next) == 0);
            Assert.assertTrue(next.equals(next));
            Assert.assertFalse(next.equals(current));
            Assert.assertFalse(current.equals(next));
            Assert.assertEquals(next.index(), current.index() + 1);
            Assert.assertEquals(next.ploidy(), ploidy);
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

    private static final int[] PLOIDY = new int[] { 1, 2, 3, 7, 10};

    @DataProvider(name="ploidyData")
    public Object[][] ploidyData() {
        final Object[][] result = new Object[PLOIDY.length][];
        for (int i = 0; i < PLOIDY.length; i++)
            result[i] = new Object[] { PLOIDY[i ]};
        return result;
    }
}
