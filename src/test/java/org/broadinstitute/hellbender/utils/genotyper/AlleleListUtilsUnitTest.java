package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Test {@link org.broadinstitute.hellbender.utils.genotyper.AlleleListUtils}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class AlleleListUtilsUnitTest {

    @Test
    public void testEmptyList(){
        final AlleleList<Allele> al = AlleleListUtils.emptyList();
        Assert.assertEquals(al.alleleCount(), 0);
        Assert.assertEquals(al.alleleIndex(null), -1);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIndexOfReferenceNull(){
        AlleleListUtils.indexOfReference(null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testasListNull(){
        AlleleListUtils.asList(null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEmptyListAdd(){
        final AlleleList<Allele> al = AlleleListUtils.emptyList();
        al.alleleAt(0);
    }

    @Test(expectedExceptions = IllegalArgumentException.class, dataProvider = "singleAlleleListData")
    public void testEqualToNull(final List<Allele> alleles1){
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles1);
        AlleleListUtils.equals(alleleList, null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class, dataProvider = "singleAlleleListData")
    public void testEqualToNull1(final List<Allele> alleles1){
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles1);
        AlleleListUtils.equals(null, alleleList);
    }

    @Test(dataProvider = "singleAlleleListData")
    public void testAsList(final List<Allele> alleles1) {
         final Allele[] uniqueAlleles = new LinkedHashSet<>(alleles1).toArray(new Allele[0]);
         final AlleleList<Allele> alleleList = new IndexedAlleleList<>(alleles1);
         final List<Allele> asList = AlleleListUtils.asList(alleleList);
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
            Assert.assertEquals(AlleleListUtils.indexOfReference(alleleList), i);
        }
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(uniqueAlleles);
        Assert.assertEquals(AlleleListUtils.indexOfReference(alleleList), -1);
    }

    @Test(dataProvider = "twoAlleleListData", dependsOnMethods={"testAsList"})
    public void testEquals(final List<Allele> alleles1, final List<Allele> alleles2) {
        final AlleleList<Allele> alleleList1 = new IndexedAlleleList<>(alleles1);
        final AlleleList<Allele> alleleList2 = new IndexedAlleleList<>(alleles2);
        Assert.assertTrue(AlleleListUtils.equals(alleleList1, alleleList1));
        Assert.assertTrue(AlleleListUtils.equals(alleleList2, alleleList2));
        Assert.assertEquals(AlleleListUtils.equals(alleleList1, alleleList2),
                Arrays.equals(AlleleListUtils.asList(alleleList1).toArray(new Allele[alleleList1.alleleCount()]),
                        AlleleListUtils.asList(alleleList2).toArray(new Allele[alleleList2.alleleCount()]))
        );
        Assert.assertEquals(AlleleListUtils.equals(alleleList1, alleleList2),
                AlleleListUtils.equals(alleleList2, alleleList1));
    }

    @Test(dataProvider = "singleAlleleListData", dependsOnMethods= "testEquals" )
    public void testSelfPermutation(final List<Allele> alleles1) {
        final AlleleList<Allele> originalAlleleList = new IndexedAlleleList<>(alleles1);
        final AlleleListPermutation<Allele> selfPermutation = AlleleListUtils.permutation(originalAlleleList,originalAlleleList);
        Assert.assertEquals(selfPermutation.fromSize(), originalAlleleList.alleleCount());
        Assert.assertEquals(selfPermutation.toSize(), originalAlleleList.alleleCount());
        Assert.assertTrue(selfPermutation.isNonPermuted());
        Assert.assertFalse(selfPermutation.isPartial());
        for (int i = 0; i < originalAlleleList.alleleCount(); i++) {
            Assert.assertEquals(selfPermutation.alleleAt(i), originalAlleleList.alleleAt(i));

            Assert.assertEquals(selfPermutation.fromIndex(i), i);
            Assert.assertEquals(selfPermutation.toIndex(i), i);
            Assert.assertEquals(selfPermutation.fromList(), selfPermutation.toList());
            AlleleListUnitTester.assertAlleleList(originalAlleleList, selfPermutation.fromList());
        }
        Assert.assertTrue(AlleleListUtils.equals(selfPermutation, originalAlleleList));
    }

    private final Random rnd = Utils.getRandomGenerator();

    @Test(dataProvider = "singleAlleleListData", dependsOnMethods = "testEquals")
    public void testSubsetPermutation(final List<Allele> alleles1) {
        final List<Allele> subsetAlleles = new ArrayList<>(alleles1.size());
        for (final Allele allele : alleles1)
            if (rnd.nextBoolean()) subsetAlleles.add(allele);
        final AlleleList<Allele> originalAlleleList = new IndexedAlleleList<>(alleles1);
        final AlleleList<Allele> targetAlleleList = new IndexedAlleleList<>(subsetAlleles);
        final AlleleListPermutation<Allele> subset = AlleleListUtils.permutation(originalAlleleList,targetAlleleList);
        if (originalAlleleList.alleleCount() == targetAlleleList.alleleCount())
            throw new SkipException("no real subset");
        Assert.assertTrue(subset.isPartial());
        Assert.assertFalse(subset.isNonPermuted());
        Assert.assertEquals(subset.fromSize(), originalAlleleList.alleleCount());
        Assert.assertEquals(subset.toSize(), targetAlleleList.alleleCount());
        AlleleListUnitTester.assertAlleleList(originalAlleleList,subset.fromList());
        AlleleListUnitTester.assertAlleleList(targetAlleleList,subset.toList());

        for (int i = 0; i < targetAlleleList.alleleCount(); i++)
            Assert.assertEquals(subset.fromIndex(i), originalAlleleList.alleleIndex(targetAlleleList.alleleAt(i)));

        for (int j = 0; j < originalAlleleList.alleleCount(); j++) {
            final Allele allele = originalAlleleList.alleleAt(j);
            Assert.assertEquals(subset.toIndex(j), targetAlleleList.alleleIndex(allele));
        }

        Assert.assertTrue(AlleleListUtils.equals(subset, targetAlleleList));
    }

    @Test(dataProvider = "singleAlleleListData", dependsOnMethods = {"testAsList","testEquals"})
    public void testShufflePermutation(final List<Allele> alleles1) {
        final AlleleList<Allele> originalAlleleList = new IndexedAlleleList<>(alleles1);
        if (originalAlleleList.alleleCount() <= 1)
            throw new SkipException("non-shuffle allele-list");

        final Allele[] targetAlleleArray = AlleleListUtils.asList(originalAlleleList).toArray(new Allele[originalAlleleList.alleleCount()]);
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

        final AlleleListPermutation<Allele> permutation = AlleleListUtils.permutation(originalAlleleList,targetAlleleList);
        Assert.assertFalse(permutation.isNonPermuted());
        AlleleListUnitTester.assertAlleleList(originalAlleleList,permutation.fromList());
        AlleleListUnitTester.assertAlleleList(targetAlleleList,permutation.toList());
        Assert.assertFalse(permutation.isPartial());
        Assert.assertEquals(permutation.fromSize(), originalAlleleList.alleleCount());
        Assert.assertEquals(permutation.toSize(), targetAlleleList.alleleCount());
        for (int i = 0; i < permutation.fromSize(); i++) {
            Assert.assertEquals(permutation.toIndex(i), targetAlleleList.alleleIndex(originalAlleleList.alleleAt(i)));
            Assert.assertEquals(permutation.fromIndex(i), originalAlleleList.alleleIndex(targetAlleleList.alleleAt(i)));
            Assert.assertEquals(permutation.fromIndex(i), fromIndex[i]);
        }
        Assert.assertTrue(AlleleListUtils.equals(permutation, targetAlleleList));

    }


    private List<Allele>[] alleleLists;

    @BeforeClass
    @SuppressWarnings("unchecked")
    public void setUp() {
        alleleLists = (List<Allele>[]) new List<?>[ALLELE_COUNT.length * MAX_ALLELE_LENGTH.length];
        int nextIndex = 0;
        for (int i = 0; i < ALLELE_COUNT.length; i++)
            for (int j = 0; j < MAX_ALLELE_LENGTH.length; j++)
                alleleLists[nextIndex++] = Arrays.asList(AlleleListUnitTester.generateRandomAlleles(ALLELE_COUNT[i], MAX_ALLELE_LENGTH[j]));
    }

    private static final int[] ALLELE_COUNT = { 0, 1, 5, 10, 20};

    private static final int[] MAX_ALLELE_LENGTH = { 1, 2, 3, 10 };

    @DataProvider(name="singleAlleleListData")
    public Object[][] singleAlleleListData() {
        final Object[][] result = new Object[alleleLists.length][];
        for (int i = 0; i < alleleLists.length; i++)
            result[i] = new Object[] { alleleLists[i]};
        return result;
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
