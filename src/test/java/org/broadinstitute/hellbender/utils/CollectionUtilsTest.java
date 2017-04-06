package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link CollectionUtils}.
 */
public class CollectionUtilsTest extends BaseTest {

    @Test(dataProvider = "heapSortingData")
    public void testHeapSorting(final List<Integer> values) {
        final List<Integer> expected = new ArrayList<>(values);
        Collections.sort(expected);
        final ArrayList<Integer> actual = CollectionUtils.heapSort(values, ArrayList::new);
        Assert.assertEquals(actual, expected);
    }

    @Test(dataProvider = "heapSortingData")
    public void testHeapSortingWithCustomComparator(final List<Integer> values) {
        // odd number then even numbers.
        // within odd numbers regular order.
        // within even number reverse order.
        final Comparator<Integer> comparator = (a, b) -> {
            if ((a & 1) == (b & 1))
                return a.compareTo(b) * (((a & 1) == 0) ? 1 : -1);
            else
                return (a & 1) == 0 ? 1 : -1;
        };
        final List<Integer> expected = new ArrayList<>(values);
        Collections.sort(expected, comparator);
        final ArrayList<Integer> actual = CollectionUtils.heapSort(values, comparator, ArrayList::new);
        Assert.assertEquals(actual, expected);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHeapSortingNullElements() {
        CollectionUtils.heapSort(Arrays.asList(1, 3, null, 5), ArrayList::new);
    }

    @DataProvider(name = "heapSortingData")
    public Object[][] heapSortingData() {
        final Random rdn = new Random(13);
        final List<Object[]> result = new ArrayList<>();
        result.add(new Object[] {Collections.emptyList()});
        result.add(new Object[] {
               IntStream.range(0, 100)
                 .map(i -> rdn.nextInt(10) * (rdn.nextBoolean() ? 0 : -1))
                 .boxed()
                 .collect(Collectors.toList())
        });
        return result.toArray(new Object[result.size()][]);
    }
}
