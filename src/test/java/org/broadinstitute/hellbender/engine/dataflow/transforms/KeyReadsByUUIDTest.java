package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestData;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.UUID;

public final class KeyReadsByUUIDTest {

    @DataProvider(name = "keyedVariantShardsReads")
    public Object[][] keyedReads(){
        DataflowTestData testData = new DataflowTestData();

        List<GATKRead> reads = testData.getReads();
        List<KV<UUID, GATKRead>> expected = Arrays.asList(
                KV.of(reads.get(0).getUUID(), reads.get(0)),
                KV.of(reads.get(1).getUUID(), reads.get(1)),
                KV.of(reads.get(2).getUUID(), reads.get(2)),
                KV.of(reads.get(3).getUUID(), reads.get(3))
        );
        return new Object[][]{
                {reads, expected},
        };
    }

    @Test(dataProvider = "keyedVariantShardsReads")
    public void selfKeyReadsTest(List<GATKRead> reads, List<KV<UUID, GATKRead>> expected) {
        Pipeline p = TestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<GATKRead> pReads = DataflowTestUtils.PCollectionCreateAndVerify(p, reads);

        PCollection<KV<UUID, GATKRead>> kReads = pReads.apply(new KeyReadsByUUID());
        DataflowAssert.that(kReads).containsInAnyOrder(expected);
        p.run();
    }
}