package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.engine.dataflow.ReadsPreprocessingPipelineTestData;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;

public final class ReferenceShardUnitTest extends BaseTest {

    @DataProvider(name = "reads")
    public Object[][] reads() {
        return new Object[][]{
                {ReadsPreprocessingPipelineTestData.makeRead(1, 300, 1, Read.class), new ReferenceShard(0,"1")},  //right in the middle of the shard
                {ReadsPreprocessingPipelineTestData.makeRead(ReferenceShard.REFERENCE_SHARD_SIZE, 10, 2, Read.class), new ReferenceShard(1, "1")}, //at  the start of a shard
                {ReadsPreprocessingPipelineTestData.makeRead(3 * ReferenceShard.REFERENCE_SHARD_SIZE - 1, 2, 3, Read.class),  new ReferenceShard(2, "1")}  //overlapping the end of a shard
        };
    }

    @Test(dataProvider = "reads")
    public void getVariantShardsFromIntervalTest(GATKRead r, ReferenceShard expectedShard) {
        ReferenceShard foundShard = ReferenceShard.getShardNumberFromInterval(r);
        Assert.assertEquals(foundShard, expectedShard);
    }

    // TODO: Add a test that verifies the coder is registered with registerGATKCoders. This will be possible with a
    // version bump of Dataflow. The new version allows users to get the coder only using the class name, not a class
    // instance.
    @Test
    public void coderTest() {
        Pipeline p = GATKTestPipeline.create();
        List<ReferenceShard> refShards = Lists.newArrayList(new ReferenceShard(0, "1"), new ReferenceShard(1, "2"), new ReferenceShard(0, "2"));
        DataflowTestUtils.pCollectionCreateAndVerify(p, refShards, SerializableCoder.of(ReferenceShard.class));
        p.run();
    }
}