package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.api.services.genomics.model.Read;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.engine.dataflow.ReadsPreprocessingPipelineTestData;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.engine.dataflow.datasources.VariantShard;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

public final class KeyReadsByOverlappingVariantShardTest extends BaseTest {

    @DataProvider(name = "keyedVariantShardsReads")
    public Object[][] keyedVariantShardsReads(){
        Object[][] data = new Object[2][];
        List<Class<?>> classes = Arrays.asList(Read.class, SAMRecord.class);
        for (int i = 0; i < classes.size(); ++i) {
            Class<?> c = classes.get(i);
            ReadsPreprocessingPipelineTestData testData = new ReadsPreprocessingPipelineTestData(c);

            List<GATKRead> reads = testData.getReads();
            List<KV<VariantShard, GATKRead>> kvVariantShardRead = testData.getKvVariantShardRead();
            data[i] = new Object[]{reads, kvVariantShardRead};
        }
        return data;
    }

    @Test(dataProvider = "keyedVariantShardsReads")
    public void keyReadsByVariantShardTest(List<GATKRead> reads, List<KV<VariantShard, GATKRead>> kvVariantShardRead) {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<GATKRead> pRead = DataflowTestUtils.pCollectionCreateAndVerify(p, reads, new GATKReadCoder());

        PCollection<KV<VariantShard, GATKRead>> kVariant = pRead.apply(new KeyReadsByOverlappingVariantShard());
        DataflowAssert.that(kVariant).containsInAnyOrder(kvVariantShardRead);
        p.run();
    }
}