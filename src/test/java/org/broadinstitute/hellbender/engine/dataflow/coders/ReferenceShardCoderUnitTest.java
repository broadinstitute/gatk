package org.broadinstitute.hellbender.engine.dataflow.coders;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.engine.ReferenceShard;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.util.List;

public class ReferenceShardCoderUnitTest extends BaseTest {
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
