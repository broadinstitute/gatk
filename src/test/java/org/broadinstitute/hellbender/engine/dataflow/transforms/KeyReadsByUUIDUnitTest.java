package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.api.services.genomics.model.Read;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.tools.ReadsPreprocessingPipelineTestData;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.UUID;

public final class KeyReadsByUUIDUnitTest extends BaseTest {

    @DataProvider(name = "keyedVariantShardsReads")
    public Object[][] keyedReads(){
        Object[][] data = new Object[2][];
        List<Class<?>> classes = Arrays.asList(Read.class, SAMRecord.class);
        for (int i = 0; i < classes.size(); ++i) {
            Class<?> c = classes.get(i);
            ReadsPreprocessingPipelineTestData testData = new ReadsPreprocessingPipelineTestData(c);

            List<GATKRead> reads = testData.getReads();
            List<KV<UUID, GATKRead>> expected = Arrays.asList(
                    KV.of(reads.get(0).getUUID(), reads.get(0)),
                    KV.of(reads.get(1).getUUID(), reads.get(1)),
                    KV.of(reads.get(2).getUUID(), reads.get(2)),
                    KV.of(reads.get(3).getUUID(), reads.get(3)),
                    KV.of(reads.get(4).getUUID(), reads.get(4))
            );
            data[i] = new Object[]{reads, expected};
        }
        return data;
    }

    @Test(dataProvider = "keyedVariantShardsReads")
    public void selfKeyReadsTest(List<GATKRead> reads, List<KV<UUID, GATKRead>> expected) {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<GATKRead> pReads = DataflowTestUtils.pCollectionCreateAndVerify(p, reads, new GATKReadCoder());

        PCollection<KV<UUID, GATKRead>> kReads = pReads.apply(new KeyReadsByUUID());
        DataflowAssert.that(kReads).containsInAnyOrder(expected);
        p.run();
    }
}