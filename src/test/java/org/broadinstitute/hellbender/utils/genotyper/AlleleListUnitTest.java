package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Test {@link org.broadinstitute.hellbender.utils.genotyper.AlleleList}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class AlleleListUnitTest {

    private final Random rnd = Utils.getRandomGenerator();

    private List<Allele>[] alleleLists;

    private static final int[] ALLELE_COUNT = { 0, 1, 5, 10, 20};
    private static final int[] MAX_ALLELE_LENGTH = { 1, 2, 3, 10 };

    @BeforeClass
    @SuppressWarnings("unchecked")
    public void setUp() {
        alleleLists = (List<Allele>[]) new List<?>[ALLELE_COUNT.length * MAX_ALLELE_LENGTH.length];
        int nextIndex = 0;
        for (int i = 0; i < ALLELE_COUNT.length; i++)
            for (int j = 0; j < MAX_ALLELE_LENGTH.length; j++)
                alleleLists[nextIndex++] = Arrays.asList(AlleleListUnitTester.generateRandomAlleles(ALLELE_COUNT[i], MAX_ALLELE_LENGTH[j]));
    }

    @BeforeMethod
    public void resetRandomGenerator(){
        Utils.resetRandomGenerator();
    }

    @Test
    public void testEmptyList(){
        final AlleleList<Allele> al = AlleleList.emptyAlleleList();
        Assert.assertEquals(al.numberOfAlleles(), 0);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEmptyListIndexOfNull(){
        final AlleleList<Allele> al = AlleleList.emptyAlleleList();
        al.indexOfAllele(null);
        Assert.assertEquals(al.indexOfAllele(null), -1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEmptyListAdd(){
        final AlleleList<Allele> al = AlleleList.emptyAlleleList();
        al.getAllele(0);
    }

    @Test
    public void testEmptyIndexedList() throws Exception {
        final AlleleList<Allele> al = new IndexedAlleleList<>();
        Assert.assertEquals(al.numberOfAlleles(), 0);
        Assert.assertEquals(al.indexOfReference(), -1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class, dataProvider = "singleAlleleListData")
    public void testEqualToNull(final List<Allele> alleles1){
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles1);
        AlleleList.equals(alleleList, null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class, dataProvider = "singleAlleleListData")
    public void testEqualToNull1(final List<Allele> alleles1){
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles1);
        AlleleList.equals(null, alleleList);
    }

    @Test(dataProvider = "singleAlleleListData")
    public void testAsList(final List<Allele> alleles1) {
         final Allele[] uniqueAlleles = new LinkedHashSet<>(alleles1).toArray(new Allele[0]);
         final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles1);
         final List<Allele> asList = alleleList.asListOfAlleles();
         final Allele[] asListArray = asList.toArray(new Allele[asList.size()]);
         Assert.assertTrue(Arrays.equals(uniqueAlleles, asListArray));
    }

    @Test(dataProvider = "singleAlleleListData")
    public void testIndexOfReference(final List<Allele> alleles1) {
        final Allele[] uniqueAlleles = new LinkedHashSet<>(alleles1).toArray(new Allele[0]);
        for (int i = 0; i < uniqueAlleles.length; i++) {
            final Allele[] actualAlleles = uniqueAlleles.clone();
            actualAlleles[i] = Allele.create(actualAlleles[i].getBases(), true);
            final AlleleList<Allele> alleleList = new IndexedAlleleList<>(actualAlleles);
            Assert.assertEquals(alleleList.indexOfReference(), i);
        }
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(uniqueAlleles);
        Assert.assertEquals(alleleList.indexOfReference(), -1);
    }

    @Test(dataProvider = "twoAlleleListData", dependsOnMethods={"testAsList"})
    public void testEquals(final List<Allele> alleles1, final List<Allele> alleles2) {
        final AlleleList<Allele> alleleList1 = new IndexedAlleleList<>(alleles1);
        final AlleleList<Allele> alleleList2 = new IndexedAlleleList<>(alleles2);
        Assert.assertTrue(AlleleList.equals(alleleList1, alleleList1));
        Assert.assertTrue(AlleleList.equals(alleleList2, alleleList2));
        Assert.assertEquals(AlleleList.equals(alleleList1, alleleList2),
                Arrays.equals(alleleList1.asListOfAlleles().toArray(new Allele[alleleList1.numberOfAlleles()]),
                        alleleList2.asListOfAlleles().toArray(new Allele[alleleList2.numberOfAlleles()]))
        );
        Assert.assertEquals(AlleleList.equals(alleleList1, alleleList2),
                AlleleList.equals(alleleList2, alleleList1));
    }

    @Test(dataProvider = "singleAlleleListData", dependsOnMethods= "testEquals" )
    public void testSelfPermutation(final List<Allele> alleles1) {
        final AlleleList<Allele> originalAlleleList = new IndexedAlleleList<>(alleles1);
        final AlleleListPermutation<Allele> selfPermutation = originalAlleleList.permutation(originalAlleleList);
        Assert.assertEquals(selfPermutation.fromSize(), originalAlleleList.numberOfAlleles());
        Assert.assertEquals(selfPermutation.toSize(), originalAlleleList.numberOfAlleles());
        Assert.assertTrue(selfPermutation.isNonPermuted());
        Assert.assertFalse(selfPermutation.isPartial());
        for (int i = 0; i < originalAlleleList.numberOfAlleles(); i++) {
            Assert.assertEquals(selfPermutation.getAllele(i), originalAlleleList.getAllele(i));

            Assert.assertEquals(selfPermutation.fromIndex(i), i);
            Assert.assertEquals(selfPermutation.toIndex(i), i);
            Assert.assertEquals(selfPermutation.fromList(), selfPermutation.toList());
            AlleleListUnitTester.assertAlleleList(originalAlleleList, selfPermutation.fromList());
        }
        Assert.assertTrue(AlleleList.equals(selfPermutation, originalAlleleList));
    }

    @Test(dataProvider = "singleAlleleListData", dependsOnMethods = "testEquals")
    public void testSubsetPermutation(final List<Allele> alleles1) {
        final List<Allele> randomSubsetAlleles = new ArrayList<>(alleles1.size());
        for (final Allele allele : alleles1) {
            if (rnd.nextBoolean()) {
                randomSubsetAlleles.add(allele);
            }
        }
        final AlleleList<Allele> originalAlleleList = new IndexedAlleleList<>(alleles1);
        final AlleleList<Allele> targetRandomAlleleList = new IndexedAlleleList<>(randomSubsetAlleles);
        final AlleleListPermutation<Allele> subset = originalAlleleList.permutation(targetRandomAlleleList);
        if (originalAlleleList.numberOfAlleles() == targetRandomAlleleList.numberOfAlleles()) {
            return;//return because this input is invalid for this test
        }
        Assert.assertTrue(subset.isPartial());
        Assert.assertFalse(subset.isNonPermuted());
        Assert.assertEquals(subset.fromSize(), originalAlleleList.numberOfAlleles());
        Assert.assertEquals(subset.toSize(), targetRandomAlleleList.numberOfAlleles());
        AlleleListUnitTester.assertAlleleList(originalAlleleList, subset.fromList());
        AlleleListUnitTester.assertAlleleList(targetRandomAlleleList, subset.toList());

        for (int i = 0; i < targetRandomAlleleList.numberOfAlleles(); i++){
            Assert.assertEquals(subset.fromIndex(i), originalAlleleList.indexOfAllele(targetRandomAlleleList.getAllele(i)));
	}

        for (int j = 0; j < originalAlleleList.numberOfAlleles(); j++) {
            final Allele allele = originalAlleleList.getAllele(j);
            Assert.assertEquals(subset.toIndex(j), targetRandomAlleleList.indexOfAllele(allele));
        }

        Assert.assertTrue(AlleleList.equals(subset, targetRandomAlleleList));
    }

    @Test(dataProvider = "singleAlleleListData", dependsOnMethods = {"testAsList","testEquals"})
    public void testShufflePermutation(final List<Allele> alleles1) {
        final AlleleList<Allele> originalAlleleList = new IndexedAlleleList<>(alleles1);
        if (originalAlleleList.numberOfAlleles() <= 1) {
            return; //return because this input is invalid for this test
        }

        final Allele[] targetAlleleArray = originalAlleleList.asListOfAlleles().toArray(new Allele[originalAlleleList.numberOfAlleles()]);
        final int[] fromIndex = new int[targetAlleleArray.length];
        for (int i = 0; i < fromIndex.length; i++)
            fromIndex[i] = i;

        for (int i = 0; i < targetAlleleArray.length - 1; i++) {
            final int swapIndex = rnd.nextInt(targetAlleleArray.length - i - 1);
            final int otherIndex = fromIndex[swapIndex + i + 1];
            final Allele other = targetAlleleArray[swapIndex + i + 1];
            fromIndex[swapIndex + i + 1] = fromIndex[i];
            fromIndex[i] = otherIndex;
            targetAlleleArray[swapIndex + i + 1] = targetAlleleArray[i];
            targetAlleleArray[i] = other;
        }
        final AlleleList<Allele> targetAlleleList = new IndexedAlleleList<>(targetAlleleArray);

        final AlleleListPermutation<Allele> permutation = originalAlleleList.permutation(targetAlleleList);
        Assert.assertFalse(permutation.isNonPermuted());
        AlleleListUnitTester.assertAlleleList(originalAlleleList, permutation.fromList());
        AlleleListUnitTester.assertAlleleList(targetAlleleList, permutation.toList());
        Assert.assertFalse(permutation.isPartial());
        Assert.assertEquals(permutation.fromSize(), originalAlleleList.numberOfAlleles());
        Assert.assertEquals(permutation.toSize(), targetAlleleList.numberOfAlleles());
        for (int i = 0; i < permutation.fromSize(); i++) {
            Assert.assertEquals(permutation.toIndex(i), targetAlleleList.indexOfAllele(originalAlleleList.getAllele(i)));
            Assert.assertEquals(permutation.fromIndex(i), originalAlleleList.indexOfAllele(targetAlleleList.getAllele(i)));
            Assert.assertEquals(permutation.fromIndex(i), fromIndex[i]);
        }
        Assert.assertTrue(AlleleList.equals(permutation, targetAlleleList));

    }

    @DataProvider(name="singleAlleleListData")
    public Iterator<Object[]> singleAlleleListData() {
        final List<Object[]> result = new ArrayList<>();
        for (int i = 0; i < alleleLists.length; i++)
            result.add(new Object[]{alleleLists[i]});
        return result.iterator();
    }

    @DataProvider(name="twoAlleleListData")
    public Object[][] twoAlleleListData() {
        final Object[][] result = new Object[alleleLists.length * alleleLists.length][];
        int index = 0;
        for (int i = 0; i < alleleLists.length; i++)
            for (int j = 0; j < alleleLists.length; j++)
                result[index++] = new Object[] { alleleLists[i], alleleLists[j]};
        return result;
    }
}