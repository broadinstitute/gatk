package org.broadinstitute.hellbender.utils.iterators;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Spliterators;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;


public class ShardedIntervalIteratorUnitTest extends GATKBaseTest {
    @Test(dataProvider = "simpleData")
    public void testSimpleCount(List<SimpleInterval> intervals, int shardSize, int expectedCount) {

        final ShardedIntervalIterator shardedIntervalIterator1 = new ShardedIntervalIterator(intervals.iterator(), shardSize);
        Assert.assertEquals(StreamSupport.stream(Spliterators.spliteratorUnknownSize(shardedIntervalIterator1, 0), false).count(), expectedCount);
    }

    @Test
    public void testNonOneValue(){
        final List<SimpleInterval> intervals = new ArrayList<>(1);
        intervals.add(new SimpleInterval("1", 500, 551));

        final int shardSizeInBases = 10;
        final ShardedIntervalIterator shardedIntervalIterator1 = new ShardedIntervalIterator(intervals.iterator(), shardSizeInBases);

        final List<SimpleInterval> newIntervals = StreamSupport.stream(Spliterators.spliteratorUnknownSize(shardedIntervalIterator1, 0), false).collect(Collectors.toList());

        Assert.assertEquals(newIntervals.size(), 6);
        Assert.assertEquals(newIntervals.get(0).size(), shardSizeInBases);
        Assert.assertEquals(newIntervals.get(0).getStart(), intervals.get(0).getStart());
        Assert.assertEquals(newIntervals.get(0).getEnd(), intervals.get(0).getStart() + shardSizeInBases - 1);

        // The last shard should only be of length 2 (base 550 & 551)
        final SimpleInterval lastInterval = newIntervals.get(newIntervals.size() - 1);
        Assert.assertEquals(lastInterval.size(), 2);
        Assert.assertEquals(lastInterval.getStart(), intervals.get(0).getEnd() - 1);
        Assert.assertEquals(lastInterval.getEnd(), intervals.get(0).getEnd());

    }



    @DataProvider(name="simpleData")
    public Object[][] getData() {
        final List<SimpleInterval> intervals = new ArrayList<>(2);
        intervals.add(new SimpleInterval("1", 100, 200));
        intervals.add(new SimpleInterval("1", 500, 550));
        return new Object[][] {
                {intervals, 1, 152},
                {intervals, 10, 11 + 6}
        };
    }
}
