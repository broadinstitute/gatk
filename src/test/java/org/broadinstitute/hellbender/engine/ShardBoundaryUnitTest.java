package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class ShardBoundaryUnitTest extends GATKBaseTest {
    private static final SimpleInterval OUTER = new SimpleInterval("1", 1, 10);
    private static final SimpleInterval INNER = new SimpleInterval("1", 3, 6);

    @Test
    public void testShardBoundaryOperations(){
        final ShardBoundary shardBoundary = new ShardBoundary(INNER, OUTER);
        Assert.assertEquals(shardBoundary.getContig(), "1");
        Assert.assertEquals(shardBoundary.getStart(), INNER.getStart());
        Assert.assertEquals(shardBoundary.getEnd(), INNER.getEnd());
        Assert.assertEquals(shardBoundary.getInterval(), INNER);
        Assert.assertEquals(shardBoundary.getPaddedInterval(), OUTER);
    }

    @DataProvider
    public Object[][] badInputs(){
        return new Object[][]{
                {null, null},
                {null, OUTER},
                {INNER, null},
                {OUTER, INNER}, //padding not within interval
        };
    }

    @Test(dataProvider = "badInputs", expectedExceptions = IllegalArgumentException.class)
    public void testShardBoundaryBadInputs(SimpleInterval interval, SimpleInterval paddingInterval){
        new ShardBoundary(interval, paddingInterval);
    }

    @DataProvider
    public Object[][] equality(){
        return new Object[][]{
                {new ShardBoundary(INNER, OUTER), new ShardBoundary(INNER, OUTER), true},
                {new ShardBoundary(INNER, OUTER), new ShardBoundary(new SimpleInterval("1", 2, 5), OUTER), false},
                {new ShardBoundary(INNER, OUTER), new ShardBoundary(INNER, INNER), false},
                {new ShardBoundary(INNER, OUTER), null, false}
        };
    }

    @Test(dataProvider = "equality")
    public void testEqualsAndHashCode(ShardBoundary first, ShardBoundary second, boolean expectedToBeEqual){
        Assert.assertEquals(first.equals(second), expectedToBeEqual);
        if (expectedToBeEqual){
            Assert.assertEquals(first.hashCode(), second.hashCode(), "hash codes must be equal if objects are equal");
        }
    }
}

