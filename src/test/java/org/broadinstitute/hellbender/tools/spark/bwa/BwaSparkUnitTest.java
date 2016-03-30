package org.broadinstitute.hellbender.tools.spark.bwa;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.util.*;

import static org.testng.Assert.assertEquals;

public final class BwaSparkUnitTest  extends BaseTest {

    @Test
    public void testConcat() {
        check(ImmutableList.of(ImmutableList.of()), ImmutableList.of());
        check(ImmutableList.of(ImmutableList.of("a")), ImmutableList.of("a"));
        check(ImmutableList.of(ImmutableList.of("a", "b")), ImmutableList.of("a", "b"));
        check(ImmutableList.of(ImmutableList.of("a"), ImmutableList.of("b")), ImmutableList.of("a", "b"));
        check(ImmutableList.of(ImmutableList.of(), ImmutableList.of("a"), ImmutableList.of("b")), ImmutableList.of("a", "b"));
        check(ImmutableList.of(ImmutableList.of("a"), ImmutableList.of(), ImmutableList.of("b")), ImmutableList.of("a", "b"));
        check(ImmutableList.of(ImmutableList.of("a"), ImmutableList.of("b"), ImmutableList.of()), ImmutableList.of("a", "b"));
        check(ImmutableList.of(ImmutableList.of("a", "b"), ImmutableList.of("c", "d")),
                ImmutableList.of("a", "b", "c", "d"));
    }

    private <T> void check(List<? extends Iterable<T>> input, List<T> expected) {
        assertEquals(Lists.newArrayList(BwaSpark.concat(input.iterator())), expected);
    }

    @Test
    public void testTransformParallel() {
        Iterator<Integer> integers = BwaSpark.transformParallel(ImmutableList.of(5, 4, 3, 2, 1).iterator(), i -> {
            try { Thread.sleep(i * 100); } catch (InterruptedException e) { }
            return i;
        }, 2);
        assertEquals(Lists.newArrayList(integers), ImmutableList.of(5, 4, 3, 2, 1));
    }
}
