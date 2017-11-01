package org.broadinstitute.hellbender.utils.collections;

import com.google.common.collect.Sets;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.LinkedHashSet;

public class NestedIntegerArrayUnitTest extends GATKBaseTest {

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testEmpty() throws Exception {
        NestedIntegerArray<String> arr= new NestedIntegerArray<>();
    }

    @Test
    public void testNew() throws Exception {
        NestedIntegerArray<String> arr= new NestedIntegerArray<>(2);
        Assert.assertEquals(arr.getDimensions(), new int[]{2});
        Assert.assertTrue(arr.getAllLeaves().isEmpty());
        Assert.assertTrue(arr.getAllValues().isEmpty());
    }

    @Test
    public void testPut() throws Exception {
        NestedIntegerArray<String> arr= new NestedIntegerArray<>(2);
        arr.put("fred", 0);
        Assert.assertEquals(arr.getDimensions(), new int[]{2});
        Assert.assertEquals(1, arr.getAllLeaves().size());
        Assert.assertEquals(Arrays.asList("fred"), arr.getAllValues());
    }

    @Test(expectedExceptions = ArrayIndexOutOfBoundsException.class)
    public void testPutBlowupTooHighValue() throws Exception {
        NestedIntegerArray<String> arr= new NestedIntegerArray<>(2);
        arr.put("fred", 3);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testPut2BlowupTooHighValue() throws Exception {
        NestedIntegerArray<String> arr= new NestedIntegerArray<>(2, 3);
        arr.put("fred", 15, 23);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testPutBlowupTooManyArgs() throws Exception {
        NestedIntegerArray<String> arr= new NestedIntegerArray<>(2);
        arr.put("fred", 0, 1);
    }

    @Test
    public void testPutget1Key() throws Exception {
        NestedIntegerArray<String> arr= new NestedIntegerArray<>(2);
        arr.put("fred", 0);
        arr.put("bozo", 1);

        Assert.assertEquals(2, arr.getAllLeaves().size());
        Assert.assertEquals(Sets.newHashSet("fred", "bozo"), new LinkedHashSet<>(arr.getAllValues()));

        Assert.assertEquals("fred", arr.get(0));
        Assert.assertEquals("fred", arr.get1Key(0));

        Assert.assertEquals("bozo", arr.get(1));
        Assert.assertEquals("bozo", arr.get1Key(1));

    }

    @Test
    public void testPutget2Keys() throws Exception {
        NestedIntegerArray<String> arr= new NestedIntegerArray<>(2, 2);
        arr.put("fred", 0, 1);
        arr.put("bozo", 1, 0);

        Assert.assertEquals(2, arr.getAllLeaves().size());
        Assert.assertEquals(Sets.newHashSet("fred", "bozo"), new LinkedHashSet<>(arr.getAllValues()));

        Assert.assertEquals("fred", arr.get(0, 1));
        Assert.assertEquals("fred", arr.get2Keys(0, 1));

        Assert.assertEquals("bozo", arr.get(1, 0));
        Assert.assertEquals("bozo", arr.get2Keys(1, 0));

        Assert.assertNull(arr.get(0, 0));
        Assert.assertNull(arr.get2Keys(0, 0));

    }

    @Test
    public void testPutget2Keysnull() throws Exception {
        NestedIntegerArray<String> arr= new NestedIntegerArray<>(2, 2);
        arr.put("fred", 0, 1);
        arr.put("bozo", 1, 0);
        arr.put(null, 1, 1);

        Assert.assertEquals(2, arr.getAllLeaves().size());
        Assert.assertEquals(Sets.newHashSet("fred", "bozo"), new LinkedHashSet<>(arr.getAllValues()));

        Assert.assertEquals("fred", arr.get(0, 1));
        Assert.assertEquals("fred", arr.get2Keys(0, 1));

        Assert.assertEquals("bozo", arr.get(1, 0));
        Assert.assertEquals("bozo", arr.get2Keys(1, 0));
    }

    @Test
    public void testPutget3Keys() throws Exception {
        NestedIntegerArray<String> arr= new NestedIntegerArray<>(2, 2, 42);
        arr.put("fred", 0, 1, 13);
        arr.put("bozo", 1, 0, 0);
        arr.put("mike", 1, 0, 17);

        Assert.assertEquals(3, arr.getAllLeaves().size());
        Assert.assertEquals(Sets.newHashSet("fred", "bozo", "mike"), new LinkedHashSet<>(arr.getAllValues()));

        Assert.assertEquals("fred", arr.get(0, 1, 13));
        Assert.assertEquals("fred", arr.get3Keys(0, 1, 13));

        Assert.assertEquals("bozo", arr.get(1, 0, 0));
        Assert.assertEquals("bozo", arr.get3Keys(1, 0, 0));

        Assert.assertEquals("mike", arr.get(1, 0, 17));
        Assert.assertEquals("mike", arr.get3Keys(1, 0, 17));

        Assert.assertNull(arr.get(0, 0, 0));
        Assert.assertNull(arr.get3Keys(0, 0, 0));

        Assert.assertNull(arr.get(1, 2, 3));
        Assert.assertNull(arr.get3Keys(1, 2, 3));
    }

    @Test
    public void testPutget4Keys() throws Exception {
        NestedIntegerArray<String> arr= new NestedIntegerArray<>(2, 2, 42, 91);
        arr.put("fred", 0, 1, 13, 41);
        arr.put("bozo", 1, 0, 0, 90);
        arr.put("mike", 1, 0, 17, 0);

        Assert.assertEquals(3, arr.getAllLeaves().size());
        Assert.assertEquals(Sets.newHashSet("fred", "bozo", "mike"), new LinkedHashSet<>(arr.getAllValues()));

        Assert.assertEquals("fred", arr.get(0, 1, 13, 41));
        Assert.assertEquals("fred", arr.get4Keys(0, 1, 13, 41));

        Assert.assertEquals("bozo", arr.get(1, 0, 0, 90));
        Assert.assertEquals("bozo", arr.get4Keys(1, 0, 0, 90));

        Assert.assertEquals("mike", arr.get(1, 0, 17, 0));
        Assert.assertEquals("mike", arr.get4Keys(1, 0, 17, 0));

        Assert.assertNull(arr.get(1, 2, 3, 4));
        Assert.assertNull(arr.get4Keys(1, 2, 3, 4));

        Assert.assertNull(arr.get(0, 0, 0, 0));
        Assert.assertNull(arr.get4Keys(0, 0, 0, 0));

    }
}

