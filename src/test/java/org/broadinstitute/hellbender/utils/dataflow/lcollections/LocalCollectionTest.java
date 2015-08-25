package org.broadinstitute.hellbender.utils.dataflow.lcollections;

import com.google.common.collect.Lists;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;

import static org.testng.Assert.*;

/**
 *
 */
public class LocalCollectionTest {

    @Test
    public void testUnion() throws Exception {
        Integer[] even = new Integer[] { 2 };
        Integer[] odd = new Integer[] { 1,3 };
        Integer[] expected =     new Integer[] { 2,1,3 };

        LocalCollection<Integer> evens = LocalCollection.of(Arrays.asList(even));
        LocalCollection<Integer> odds = LocalCollection.of(Arrays.asList(odd));
        LocalCollection<Integer> together = LocalCollection.union(evens, odds);

        assertEquals(expected, together.iterable());

    }

    @Test
    public void testPartitionBy() throws Exception {
        Integer[] nums =     new Integer[] { 1,2,3 };
        Integer[] even = new Integer[] { 2 };
        Integer[] odd = new Integer[] { 1,3 };
        LocalCollection<Integer> col = LocalCollection.of(Arrays.asList(nums));
        ArrayList<LocalCollection<Integer>> parity = col.partitionBy((i) -> (i % 2));
        Assert.assertEquals(parity.size(), 2);
        assertEquals(even, parity.get(0).iterable());
        assertEquals(odd, parity.get(1).iterable());
    }

    @Test
    public void testPartitionDisjoint() throws Exception {
        Integer[] nums =     new Integer[] { 1,2,3 };
        Integer[] even = new Integer[] { 2 };
        Integer[] odd = new Integer[] { 1,3 };
        LocalCollection<Integer> col = LocalCollection.of(Arrays.asList(nums));
        ArrayList<LocalCollection<Integer>> parity = col.partitionBy((i) -> (i % 2==0?2:4));
        Assert.assertEquals(parity.size(), 5);
        assertEquals(even, parity.get(2).iterable());
        assertEquals(odd, parity.get(4).iterable());
    }

    @Test
    public void testPartitionMinSize() throws Exception {
        Integer[] nums =     new Integer[] { 1,2,3 };
        Integer[] even = new Integer[] { 2 };
        Integer[] odd = new Integer[] { 1,3 };
        LocalCollection<Integer> col = LocalCollection.of(Arrays.asList(nums));
        ArrayList<LocalCollection<Integer>> parity = col.partitionBy((i) -> (i % 2), 10);
        Assert.assertEquals(parity.size(), 10);
        assertEquals(even, parity.get(0).iterable());
        assertEquals(odd, parity.get(1).iterable());
    }

    @Test
    public void testGroupBy() throws Exception {
        Integer[] nums =     new Integer[] { 1,2,3,4,5 };
        Integer[] even = new Integer[] { 2 };
        Integer[] odd = new Integer[] { 1,3 };
        LocalCollection<Integer> evens = LocalCollection.of(Arrays.asList(even));
        LocalCollection<Integer> odds = LocalCollection.of(Arrays.asList(odd));

        LocalGroupedCollection<Integer> col =
            LocalCollection.of(Arrays.asList(nums))
            .groupBy((i) -> "" + (i % 2));
        assertEquals(even, col.get("0").iterable());
        assertEquals(odd, col.get("1").iterable());
    }

    @Test
    public void testMap() throws Exception {
        Integer[] nums =     new Integer[] { 1,2,3 };
        Integer[] expected = new Integer[] { 2,3,4 };
        LocalCollection<Integer> col = LocalCollection.of(Arrays.asList(nums));
        col = col.map((i)->i+1);
        assertEquals(expected, col.iterable());
    }

    @Test
    public void testTransform() throws Exception {

    }

    public static <T> void assertEquals(T[] array, Iterable<T> iterable) {
        ArrayList<T> actual = Lists.newArrayList(iterable);
        Assert.assertEquals(array.length, actual.size());
        for (int i=0; i<array.length; i++) {
            Assert.assertEquals( array[i], actual.get(i));
        }
    }
}