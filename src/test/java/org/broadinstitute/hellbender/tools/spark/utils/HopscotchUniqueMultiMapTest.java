package org.broadinstitute.hellbender.tools.spark.utils;

import org.broadinstitute.hellbender.tools.spark.utils.HopscotchMapTest.IntPair;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class HopscotchUniqueMultiMapTest extends GATKBaseTest {
    // a HopscotchUniqueMultiMap is just like a HopscotchMultiMap, which is tested separately, except that it
    // cannot contain multiple equivalent entries.
    // multiple entries for a given key, yes.  multiple entries equal to each other, no.
    // so we just test that feature

    @Test
    void noDupsTest() {
        final HopscotchUniqueMultiMap<Integer, Integer, IntPair> multiMap = new HopscotchUniqueMultiMap<>();
        multiMap.add(new IntPair(1, 0));
        Assert.assertFalse(multiMap.add(new IntPair(1, 0)));
        Assert.assertEquals(multiMap.size(), 1);
    }
}
