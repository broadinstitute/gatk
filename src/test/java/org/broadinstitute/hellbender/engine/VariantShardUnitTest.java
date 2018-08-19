package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.testutils.ReadsPreprocessingPipelineTestData;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

public final class VariantShardUnitTest extends GATKBaseTest {

    @DataProvider(name = "variantShards")
    public Object[][] reads() {
        return new Object[][]{
                {ReadsPreprocessingPipelineTestData.makeRead("1", 1, 300, 1, SAMRecord.class),
                        Arrays.asList(new VariantShard(0, "1"))},  //right in the middle of the shard
                {ReadsPreprocessingPipelineTestData.makeRead("1", VariantShard.VARIANT_SHARDSIZE, 10, 2, SAMRecord.class),
                        Arrays.asList(new VariantShard(1, "1"))}, //at  the start of a shard
                {ReadsPreprocessingPipelineTestData.makeRead("1", 3 * VariantShard.VARIANT_SHARDSIZE - 1, 2, 3, SAMRecord.class),
                        Arrays.asList(new VariantShard(2, "1"), new VariantShard(3, "1"))}  // overlapping the end of a shard
        };
    }

    @Test(dataProvider = "variantShards")
    public void getVariantShardsFromIntervalTest(GATKRead read, Iterable<VariantShard> expectedShards) {
            List<VariantShard> foundShards = VariantShard.getVariantShardsFromInterval(read);
            Assert.assertEquals(foundShards, expectedShards);
    }
}