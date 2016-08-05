package org.broadinstitute.hellbender.engine.spark;

import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.SparkTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class SparkReadShardUnitTest extends BaseTest {
    private static final SimpleInterval OUTER = new SimpleInterval("1", 1, 10);
    private static final SimpleInterval INNER = new SimpleInterval("1", 5, 7);
    private static final ShardBoundary BOUND = new ShardBoundary(INNER, OUTER);

    @Test
    public void testNullValues(){
        Assert.assertThrows(IllegalArgumentException.class, () -> new SparkReadShard(null, Collections.emptyList()));
        Assert.assertThrows(IllegalArgumentException.class, () -> new SparkReadShard(BOUND, null));
    }

    @Test
    public void testSparkReadShardBasicOperations(){
        final List<GATKRead> readList = Arrays.asList(ArtificialReadUtils.createArtificialRead("50M"));
        final SparkReadShard shard = new SparkReadShard(BOUND, readList);
        Assert.assertEquals(shard.getInterval(), INNER);
        Assert.assertEquals(shard.getPaddedInterval(), OUTER);
        Assert.assertEquals(shard.getContig(), "1");
        Assert.assertEquals(shard.getStart(), INNER.getStart());
        Assert.assertEquals(shard.getEnd(), INNER.getEnd());
        Assert.assertEquals(shard.iterator().next(), readList.get(0));
    }

    @Test
    public void testKryoCanSerialize(){
        final List<GATKRead> readList = Arrays.asList(ArtificialReadUtils.createArtificialRead("50M"));
        final SparkReadShard initialShard = new SparkReadShard(BOUND, readList);
        final SparkReadShard shard = SparkTestUtils.roundTripInKryo(initialShard, SparkReadShard.class, SparkContextFactory.getTestSparkContext().getConf());
        Assert.assertEquals(shard, initialShard);
    }
}

