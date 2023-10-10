package org.broadinstitute.hellbender.utils;

import com.google.common.collect.Streams;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.List;

public class SmallBitSetUnitTest {
    private static final SmallBitSet EMPTY = new SmallBitSet();
    private static final int ELEMENT_1 = 5;
    private static final int ELEMENT_2 = 13;
    private static final int ELEMENT_3 = 7;
    private static final int DIFFERENT_ELEMENT = 20;
    private static final SmallBitSet ONE_ELEMENT = new SmallBitSet(ELEMENT_1);
    private static final SmallBitSet TWO_ELEMENTS = new SmallBitSet(ELEMENT_1, ELEMENT_2);
    private static final SmallBitSet THREE_ELEMENTS = new SmallBitSet(ELEMENT_1, ELEMENT_2, ELEMENT_3);

    @Test
    public void testConstructors() {
        for (int i = 0; i < SmallBitSet.MAX_ELEMENTS; i++) {
            Assert.assertFalse(EMPTY.get(i));
            Assert.assertEquals(i == ELEMENT_1, ONE_ELEMENT.get(i));
            Assert.assertEquals(i == ELEMENT_1 || i == ELEMENT_2, TWO_ELEMENTS.get(i));
            Assert.assertEquals(i == ELEMENT_1 || i == ELEMENT_2 || i == ELEMENT_3, THREE_ELEMENTS.get(i));
        }
    }

    @Test
    public void testIntersection() {
        Assert.assertEquals((new SmallBitSet(1,3,5)).intersection(new SmallBitSet(7,9,11)), EMPTY);

        Assert.assertEquals((new SmallBitSet(1,3,5)).intersection(new SmallBitSet(5,7,9)),
                (new SmallBitSet(5)));

        Assert.assertEquals((new SmallBitSet(List.of(0, 20, 5, 15, 10)))
                .intersection(new SmallBitSet(List.of(0, 7, 21, 14, 20))), new SmallBitSet(0, 20));
    }

    @Test
    public void testUnion() {
        Assert.assertEquals((new SmallBitSet(1,3,5)).union(new SmallBitSet(7,9,11)),
                new SmallBitSet(List.of(1,3,5,7,9,11)));

        Assert.assertEquals((new SmallBitSet(1,3,5)).union(new SmallBitSet(5,7,9)),
                new SmallBitSet(List.of(1,3,5,7,9)));
    }

    @Test
    public void testStream() {
        Assert.assertEquals((new SmallBitSet(1,3,5)).stream(10).toArray(),
                new int[] {1,3,5});

        Assert.assertEquals((new SmallBitSet(1,5,11)).stream(10).toArray(),
                new int[] {1, 5});
    }

    @Test
    public void testContains() {
        Assert.assertTrue(THREE_ELEMENTS.contains(TWO_ELEMENTS));
        Assert.assertTrue(THREE_ELEMENTS.contains(ONE_ELEMENT));
        Assert.assertTrue(TWO_ELEMENTS.contains(ONE_ELEMENT));
        Assert.assertTrue(THREE_ELEMENTS.contains(THREE_ELEMENTS));
        Assert.assertTrue(TWO_ELEMENTS.contains(TWO_ELEMENTS));
        Assert.assertTrue(ONE_ELEMENT.contains(ONE_ELEMENT));
        Assert.assertFalse(TWO_ELEMENTS.contains(THREE_ELEMENTS));
        Assert.assertFalse(ONE_ELEMENT.contains(THREE_ELEMENTS));
        Assert.assertFalse(ONE_ELEMENT.contains(TWO_ELEMENTS));

        Assert.assertFalse(THREE_ELEMENTS.contains(new SmallBitSet(ELEMENT_1, ELEMENT_2, DIFFERENT_ELEMENT)));
    }

    @Test
    public void testAdd() {
        final SmallBitSet copy = EMPTY.copy();
        copy.add(ELEMENT_1);
        Assert.assertEquals(copy, ONE_ELEMENT);
        copy.add(ELEMENT_2);
        Assert.assertEquals(copy, TWO_ELEMENTS);
        copy.add(ELEMENT_3);
        Assert.assertEquals(copy, THREE_ELEMENTS);
    }

    @Test
    public void testRemove() {
        final SmallBitSet copy = THREE_ELEMENTS.copy();
        copy.remove(ELEMENT_3);
        Assert.assertEquals(copy, TWO_ELEMENTS);
        copy.remove(ELEMENT_2);
        Assert.assertEquals(copy, ONE_ELEMENT);
        copy.remove(ELEMENT_1);
        Assert.assertEquals(copy, EMPTY);
    }

    @Test
    public void testFlip() {
        final SmallBitSet bitset = new SmallBitSet(List.of(5, 10, 15, 20));
        bitset.flip(10);
        Assert.assertEquals(bitset, new SmallBitSet(5, 15, 20));
        bitset.flip(7);
        Assert.assertEquals(bitset, new SmallBitSet(List.of(5, 7, 15, 20)));
    }

    @Test
    public void testGet() {
        List.of(ELEMENT_1, ELEMENT_2, ELEMENT_3).forEach(el -> Assert.assertTrue(THREE_ELEMENTS.get(el)));
        Assert.assertFalse(THREE_ELEMENTS.get(DIFFERENT_ELEMENT));
    }
}