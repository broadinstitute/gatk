package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.RandomDNA;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.SkipException;

import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

/**
 * Helper class for those unit-test classes that test on implementations of SampleList.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class AlleleListUnitTester {

    private static final Random rnd = Utils.getRandomGenerator();
    private static final RandomDNA rndDNA = new RandomDNA(rnd);

    /**
     * Test that the contents of an allele-list are the ones expected.
     * <p/>
     * <p>
     * This method perform various consistency check involving all the {@link org.broadinstitute.hellbender.utils.genotyper.AlleleList} interface methods.
     * Therefore calling this method is equivalent to a thorough check of the {@link org.broadinstitute.hellbender.utils.genotyper.AlleleList} aspect of
     * the {@code actual} argument.
     * </p>
     *
     * @param actual   the sample-list to assess.
     * @param expected the expected sample-list.
     * @throws IllegalArgumentException if {@code expected} is {@code null} or contains
     *                                  {@code null}s which is an indication of an bug in the testing code.
     * @throws RuntimeException         if there is some testing assertion exception which
     *                                  is an indication of an actual bug the code that is been tested.
     */
    @SuppressWarnings("unchecked")
    public static <A extends Allele> void assertAlleleList(final AlleleList<A> actual, final List<A> expected) {
        Utils.nonNull(expected, "the expected list cannot be null");
        final Set<A> expectedAlleleSet = new LinkedHashSet<>(expected.size());
        Assert.assertNotNull(actual);
        Assert.assertEquals(actual.numberOfAlleles(), expected.size());
        for (int i = 0; i < expected.size(); i++) {
            final A expectedAllele = expected.get(i);
            Utils.nonNull(expectedAllele, "the expected sample cannot be null");
            Utils.validateArg(!expectedAllele.equals(NEVER_USE_ALLELE), "you cannot use the forbidden sample name");
            Utils.validateArg(!expectedAlleleSet.contains(expected.get(i)), "repeated allele in the expected list, this is a test bug");
            final A actualAllele = actual.getAllele(i);
            Assert.assertNotNull(actualAllele, "allele cannot be null");
            Assert.assertFalse(expectedAlleleSet.contains(actualAllele), "repeated allele: " + actualAllele);
            Assert.assertEquals(actualAllele, expectedAllele, "wrong allele order; index = " + i);
            Assert.assertEquals(actual.indexOfAllele(actualAllele), i, "allele index mismatch");
            expectedAlleleSet.add(actualAllele);
        }

        Assert.assertEquals(actual.indexOfAllele((A) NEVER_USE_ALLELE), -1);
    }

    /**
     * Save to assume that this allele will never be used.
     */
    private static final Allele NEVER_USE_ALLELE = Allele.create(new String("ACTGACTGACTGACTGACTGACTGACTGACTGGTCAGTCAGTCAGTCAGTCAGTCA").getBytes(), false);

    /**
     * Generate testing alleles. Guarantees that alleles are unique.
     *
     * <p>
     *     Basically all are random alleles given the maximum allele length.
     * </p>
     *
     * <p>
     *     So with a low max-allele-length and high allele-count you can force repeats.
     * </p>
     *
     * @param alleleCount number of alleles to generate.
     * @param maxAlleleLength the maximum length of the allele in bases.
     *
     * @throws RuntimeException if {@code alleleCount} is negative or {@code maxAlleleLength} is less than 1.
     * @return never {@code null}.
     */
    public static Allele[] generateRandomUniqueAlleles(final int alleleCount, final int maxAlleleLength) {
        Utils.validateArg(maxAlleleLength > 0, "the max allele length cannot be less than 1");
        final Set<Allele> set = new LinkedHashSet<>(alleleCount);
        while(set.size() < alleleCount){
            final int alleleLength = rnd.nextInt(maxAlleleLength) + 1;
            final Allele result = Allele.create(rndDNA.nextBases(alleleLength));
            set.add(result);
        }
        return set.toArray(new Allele[set.size()]);
    }

    /**
     * Generate testing alleles.
     *
     * <p>
     *     Basically all are random alleles given the maximum allele length.
     * </p>
     *
     * <p>
     *     So with a low max-allele-length and high allele-count you can force repeats.
     * </p>
     *
     * @param alleleCount number of alleles to generate.
     * @param maxAlleleLength the maximum length of the allele in bases.
     *
     * @throws RuntimeException if {@code alleleCount} is negative or {@code maxAlleleLength} is less than 1.
     * @return never {@code null}.
     */
    public static Allele[] generateRandomAlleles(final int alleleCount, final int maxAlleleLength) {
        if (maxAlleleLength < 1)
            throw new IllegalArgumentException("the max allele length cannot be less than 1");
        final Allele[] result = new Allele[alleleCount];
        for (int i = 0; i < alleleCount; i++) {
            final int alleleLength = rnd.nextInt(maxAlleleLength) + 1;
            result[i] = Allele.create(rndDNA.nextBases(alleleLength));
        }
        return result;
    }

    /**
     * Generate testing alleles.
     *
     * <p>
     *     Basically all are random alleles given the maximum allele length.
     * </p>
     *
     * <p>
     *     So with a low max-allele-length and high allele-count you can force repeats.
     * </p>
     *
     * @param alleleCount number of alleles to generate.
     * @param maxAlleleLength the maximum length of the allele in bases.
     * @param skipIfRepeats throw an test-skip exception {@link SkipException} if the resulting allele-list
     *                     has repeats, thus is size is less than {@code alleleCount}
     *
     * @throws RuntimeException if {@code alleleCount} is negative or {@code maxAlleleLength} is less than 1.
     * @return never {@code null}.
     */
    public static AlleleList<Allele> alleleList(final int alleleCount, final int maxAlleleLength, final boolean skipIfRepeats) {
        final Allele[] alleles = AlleleListUnitTester.generateRandomUniqueAlleles(alleleCount, maxAlleleLength);
        if (alleleCount > 0)
            alleles[0] = Allele.create(alleles[0].getBases(), true);
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles);
        if (skipIfRepeats && alleleList.numberOfAlleles() != alleles.length)
            throw new SkipException("repeated alleles, should be infrequent");
        return alleleList;
    }
}