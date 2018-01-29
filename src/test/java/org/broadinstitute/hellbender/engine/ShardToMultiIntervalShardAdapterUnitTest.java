package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.engine.ShardToMultiIntervalShardAdapter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class ShardToMultiIntervalShardAdapterUnitTest extends GATKBaseTest {
    private static final SimpleInterval INTERVAL = new SimpleInterval("1", 10, 15);
    private static final SimpleInterval PADDED_INTERVAL = new SimpleInterval("1", 5, 22);
    private static final List<Integer> INTEGERS = Arrays.asList(1, 2, 3);
    private static final ShardToMultiIntervalShardAdapter<Integer> shard = new ShardToMultiIntervalShardAdapter<>(new Shard<Integer>(){

        @Override
        public Iterator<Integer> iterator() {
            return INTEGERS.iterator();
        }

        @Override
        public SimpleInterval getInterval() {
           return INTERVAL;
        }

        @Override
        public SimpleInterval getPaddedInterval() {
            return PADDED_INTERVAL;
        }
    });

    @Test
    public void testGetIntervals() {
        Assert.assertEquals(shard.getIntervals().get(0), INTERVAL);
        Assert.assertEquals(shard.getIntervals().size(), 1);
    }

    @Test
    public void testGetPaddedIntervals() {
        Assert.assertEquals(shard.getPaddedIntervals().get(0), PADDED_INTERVAL);
        Assert.assertEquals(shard.getPaddedIntervals().size(), 1);
    }

    @Test
    public void testIterator() {
        Assert.assertEquals(shard.iterator(), INTEGERS.iterator());
    }

    @Test
    public void testGetInterval() {
        Assert.assertEquals(shard.getInterval(), INTERVAL);
    }

    @Test
    public void testGetPaddedInterval() {
        Assert.assertEquals(shard.getPaddedInterval(), PADDED_INTERVAL);
    }
}