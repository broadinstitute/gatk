package org.broadinstitute.hellbender.utils.collections;

import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Tests the working of {@link IndexedSet}
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class IndexedSetUnitTest {

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInitError(){
        final Collection<String> c = null;
        final IndexedSet<String> is = new IndexedSet<>(c);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEqualsError(){
        final String[] c = {"a", "b"};
        final IndexedSet<String> is = new IndexedSet<>(c);
        is.equals(null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInitError_nullColl(){
        final Collection<String> c = Arrays.asList("a", null);
        final IndexedSet<String> is = new IndexedSet<>(c);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInitError2(){
        final String[] c = null;
        final IndexedSet<String> is = new IndexedSet<>(c);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInitError3(){
        final String[] c = {"a", null};
        final IndexedSet<String> is = new IndexedSet<>(c);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInitError4(){
        final String[] c = {"a", "b"};
        final IndexedSet<String> is = new IndexedSet<>(c);
        is.add(null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInitError5(){
        final String[] c = {"a", "b"};
        final IndexedSet<String> is = new IndexedSet<>(c);
        is.indexOf(null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testInitError6(){
        final String[] c = {"a", "b"};
        final IndexedSet<String> is = new IndexedSet<>(c);
        is.remove(null);
    }

    @Test
    public void testRemove(){
        final String[] c = {"a", "b"};
        final IndexedSet<String> is = new IndexedSet<>(c);
        final boolean cC = is.remove("c");
        Assert.assertEquals(cC, false);

        final boolean cb = is.remove("b");
        Assert.assertEquals(cb, true);
    }

    @Test
    public void testEquals(){
        final String[] c = {"a", "b"};
        final String[] c1 = {"a", "B"};
        final String[] c3 = {"a", "b", "c"};
        final IndexedSet<String> is1 = new IndexedSet<>(c);
        final IndexedSet<String> is2 = new IndexedSet<>(c);
        final IndexedSet<String> is3 = new IndexedSet<>(c1);
        final IndexedSet<String> is4 = new IndexedSet<>(c3);
        Assert.assertEquals(is1, is2);
        Assert.assertEquals(is1, is1);

        Assert.assertTrue(is1.equals(is1));
        Assert.assertTrue(is1.equals((Object)is1));

        final String nullString = null;
        Assert.assertFalse(is1.equals(nullString));
        Assert.assertFalse(is1.equals("a,b"));

        Assert.assertNotEquals(is1, is3);
        Assert.assertNotEquals(is2, is3);
        Assert.assertNotEquals(is1, is4);
        Assert.assertNotEquals(is2, is4);
        Assert.assertNotEquals(is3, is4);

        Assert.assertEquals(is1.hashCode(), is2.hashCode());
        Assert.assertNotEquals(is1.hashCode(), is3.hashCode());
        Assert.assertNotEquals(is2.hashCode(), is3.hashCode());
    }

    @Test(dataProvider = "initialCapacityElementCountMaxElementData")
    public void testCompositionBySingleElementAddition(final int initialCapacity,
                                                       final int elementCount, final int maxElement) {
        final Random rnd = Utils.getRandomGenerator();
        final IndexedSet<Integer> subject = new IndexedSet<>(initialCapacity);

        final Set<Integer> elementSet = new LinkedHashSet<>();

        for (int i = 0; i < elementCount; i++) {
            final int nextElement = rnd.nextInt(maxElement + 1);
            final boolean isNewElement = ! elementSet.contains(nextElement);
            Assert.assertEquals(subject.add(nextElement), elementSet.add(nextElement));
            Assert.assertEquals(subject.size(), elementSet.size());
            if (isNewElement)
                Assert.assertEquals(subject.indexOf(nextElement), elementSet.size() - 1);
        }
        assertEquals(subject, elementSet);
    }

    @Test(dataProvider = "initialCapacityElementCountMaxElementData")
    public void testCompositionByCollectionAddition(final int initialCapacity,
                                                    final int elementCount, final int maxElement) {
        final IndexedSet<Integer> subject = new IndexedSet<>(initialCapacity);
        final List<Integer> elementList = generateElementCollection(elementCount,maxElement);


        Assert.assertEquals(subject.addAll(elementList), !elementList.isEmpty());

        final Set<Integer> elementSet = new LinkedHashSet<>(elementCount);
        elementSet.addAll(elementList);

        assertEquals(subject,elementSet);
    }

    @Test(dataProvider = "elementCountMaxElementData")
    public void testCompositionByCollectionConstructor(final int elementCount, final int maxElement) {
        final List<Integer> elementList = generateElementCollection(elementCount, maxElement);

        final IndexedSet<Integer> subject = new IndexedSet<>(elementList);

        final Set<Integer> elementSet = new LinkedHashSet<>(elementList);
        assertEquals(subject,elementSet);
        Assert.assertFalse(subject.addAll(elementList));
    }

    private List<Integer> generateElementCollection(final int elementCount, final int maxElement) {
        final Random rnd = Utils.getRandomGenerator();

        final List<Integer> elementList = new ArrayList<>(elementCount);
        for (int i = 0; i < elementCount; i++)
            elementList.add(rnd.nextInt(maxElement + 1));
        return elementList;
    }

    @Test(dataProvider = "elementCountMaxElementData",
          dependsOnMethods = {"testCompositionByCollectionConstructor"})
    public void testLookupByIndex(final int elementCount, final int maxElement) {
        final List<Integer> elementList = generateElementCollection(elementCount, maxElement);
        final IndexedSet<Integer> subject = new IndexedSet<>(elementList);
        final Set<Integer> elementSet = new LinkedHashSet<>(elementList);
        final Integer[] elementArray = elementSet.toArray(new Integer[elementSet.size()]);

        final List<Integer> subjectList = subject.asList();
        for (int i = 0; i < subject.size(); i++) {
            final int element = elementArray[i];
            final int subjectElement = subject.get(i);
            final int subjectListElement = subjectList.get(i);
            Assert.assertEquals(subjectElement, element);
            Assert.assertEquals(subjectListElement, element);
        }
    }

    @Test(dataProvider = "elementCountMaxElementData",
            dependsOnMethods = {"testCompositionByCollectionConstructor"})
    public void testIndexOf(final int elementCount, final int maxElement) {
        final List<Integer> elementList = generateElementCollection(elementCount, maxElement);
        final IndexedSet<Integer> subject = new IndexedSet<>(elementList);
        final Set<Integer> elementSet = new LinkedHashSet<>(elementList);
        final Integer[] elementArray = elementSet.toArray(new Integer[elementSet.size()]);

        final List<Integer> subjectList = subject.asList();
        for (int i = 0; i < subject.size(); i++) {
            final int element = elementArray[i];
            final int listElement = subjectList.get(i);
            final int subjectIndex = subject.indexOf(element);
            Assert.assertEquals(listElement, element);
            Assert.assertEquals(subjectIndex, i);
            Assert.assertEquals(subject.indexOf(-element - 1), -1);
        }
    }

    @Test(dataProvider = "elementCountMaxElementData",
            dependsOnMethods = {"testCompositionByCollectionConstructor","testIndexOf"})
    public void testRemoveHalf(final int elementCount, final int maxElement) {
        final List<Integer> elementList = generateElementCollection(elementCount, maxElement);
        final IndexedSet<Integer> subject = new IndexedSet<>(elementList);
        final Set<Integer> elementSet = new LinkedHashSet<>(elementList);
        final int removeCount = (subject.size() + 1) / 2;
        final Random rnd = Utils.getRandomGenerator();
        for (int i = 0; i < removeCount; i++) {
            final int removeIndex = rnd.nextInt(subject.size());
            final int removeElement = subject.get(removeIndex);
            subject.remove(removeElement);
            elementSet.remove(removeElement);
        }

        assertEquals(subject,elementSet);
    }

    @Test(dataProvider = "elementCountMaxElementData",
            dependsOnMethods = {"testCompositionByCollectionConstructor","testIndexOf"})
    public void testRemoveAll(final int elementCount, final int maxElement) {
        final List<Integer> elementList = generateElementCollection(elementCount, maxElement);
        final IndexedSet<Integer> subject = new IndexedSet<>(elementList);
        final Set<Integer> elementSet = new LinkedHashSet<>(elementList);
        final int removeCount = subject.size();
        final Random rnd = Utils.getRandomGenerator();
        for (int i = 0; i < removeCount; i++) {
            final int removeIndex = rnd.nextInt(subject.size());
            final int removeElement = subject.get(removeIndex);
            subject.remove(removeElement);
            elementSet.remove(removeElement);
        }

        assertEquals(subject,elementSet);
    }

    @Test(dataProvider = "elementCountMaxElementData",
            dependsOnMethods = {"testCompositionByCollectionConstructor"})
    public void testClear(final int elementCount, final int maxElement) {
        final List<Integer> elementList = generateElementCollection(elementCount, maxElement);
        final IndexedSet<Integer> subject = new IndexedSet<>(elementList);
        final Set<Integer> elementSet = new LinkedHashSet<>(elementList);
        subject.clear();
        elementSet.clear();

        assertEquals(subject, elementSet);
    }

    @Test(dataProvider = "elementCountMaxElementData",
            dependsOnMethods = {"testCompositionByCollectionConstructor","testIndexOf"})
    public void testRemoveAndAdd(final int elementCount, final int maxElement) {
        final List<Integer> elementList = generateElementCollection(elementCount, maxElement);
        final IndexedSet<Integer> subject = new IndexedSet<>(elementList);
        final Set<Integer> elementSet = new LinkedHashSet<>(elementList);
        final int removeCount = subject.size();
        final Random rnd = Utils.getRandomGenerator();
        for (int i = 0; i < removeCount; i++) {
            final int removeIndex = rnd.nextInt(subject.size());
            final int removeElement = subject.get(removeIndex);
            subject.remove(removeElement);
            elementSet.remove(removeElement);
        }
        subject.addAll(elementList);
        elementSet.addAll(elementList);

        assertEquals(subject, elementSet);
    }

    private final int[] INITIAL_CAPACITY = { 0, 10, 100 };

    private final int[] ELEMENT_COUNT = { 0, 1, 10, 100 , 1000 };

    private final int[] MAX_ELEMENT = { 0, 1, 5, 10, 50, 100, 500 };

    @DataProvider(name="initialCapacityElementCountMaxElementData")
    public Object[][] initialCapacityElementCountMaxElementData() {
        final Object[][] result = new Object[INITIAL_CAPACITY.length * ELEMENT_COUNT.length * MAX_ELEMENT.length][];

        int nextIndex = 0;

        for (int i = 0; i < INITIAL_CAPACITY.length; i++)
            for (int j = 0; j < ELEMENT_COUNT.length; j++)
                for (int k = 0; k < MAX_ELEMENT.length; k++)
                    result[nextIndex++] = new Object[] { INITIAL_CAPACITY[i], ELEMENT_COUNT[j], MAX_ELEMENT[k]};

        return result;
    }

    @DataProvider(name="elementCountMaxElementData")
    public Object[][] elementCountMaxElementData() {
        final Object[][] result = new Object[ELEMENT_COUNT.length * MAX_ELEMENT.length][];

        int nextIndex = 0;

            for (int j = 0; j < ELEMENT_COUNT.length; j++)
                for (int k = 0; k < MAX_ELEMENT.length; k++)
                    result[nextIndex++] = new Object[] { ELEMENT_COUNT[j], MAX_ELEMENT[k]};

        return result;
    }

    /**
     * Asserts that an indexed-set is equivalent to a insertion-sorted set provided.
     * @param subject  the indexed-set to test.
     * @param elementSet the insertion-sorted set.
     */
    private void assertEquals(final IndexedSet<Integer> subject, final Set<Integer> elementSet) {
        Assert.assertEquals(subject.size(), elementSet.size());
        final List<Integer> subjectList = subject.asList();
        Assert.assertEquals(subjectList.size(), elementSet.size());
        final Iterator<Integer> subjectIterator = subject.iterator();
        final Iterator<Integer> elementSetIterator = subject.iterator();

        final ListIterator<Integer> subjectListIterator = subjectList.listIterator();

        while (subjectIterator.hasNext()) {
            Assert.assertTrue(elementSetIterator.hasNext(), "less elements in indexed-set than in the equivalent hash-set");
            Assert.assertTrue(subjectListIterator.hasNext());

            final Integer nextElement;
            Assert.assertEquals(nextElement = subjectIterator.next(), elementSetIterator.next(), "elements in indexed-set do not follow the same order as equivalent linked hash-set's");
            Assert.assertEquals(subjectListIterator.next(), nextElement);
            Assert.assertEquals(subject.indexOf(nextElement), subjectListIterator.previousIndex());
        }
        Assert.assertFalse(elementSetIterator.hasNext());
        Assert.assertFalse(subjectListIterator.hasNext());
    }
}
