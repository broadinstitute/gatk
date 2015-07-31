package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.IterableCoder;
import com.google.cloud.dataflow.sdk.coders.KvCoder;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.engine.dataflow.ReadsPreprocessingPipelineTestData;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceShard;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

public final class KeyReadsByRefShardUnitTest extends BaseTest {

    @DataProvider(name = "refShards")
    public Object[][] refShards(){
        Object[][] data = new Object[2][];
        List<Class<?>> classes = Arrays.asList(Read.class, SAMRecord.class);
        for (int i = 0; i < classes.size(); ++i) {
            Class<?> c = classes.get(i);
            ReadsPreprocessingPipelineTestData testData = new ReadsPreprocessingPipelineTestData(c);

            List<GATKRead> inputs = testData.getReads();
            List<KV<ReferenceShard, Iterable<GATKRead>>> kvs = testData.getKvRefShardiReads();
            data[i] = new Object[]{inputs, kvs};

        }
        return data;
    }

    @Test(dataProvider = "refShards")
    public void groupReadsForRefTest(List<GATKRead> reads, List<KV<ReferenceShard, Iterable<GATKRead>>> expectedResult) {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<GATKRead> pReads = DataflowTestUtils.pCollectionCreateAndVerify(p, reads, new GATKReadCoder());
        PCollection<KV<ReferenceShard, Iterable<GATKRead>>> grouped = pReads.apply(new KeyReadsByRefShard());

        PCollection<KV<ReferenceShard, Iterable<GATKRead>>> pFinalExpected = p.apply(Create.of(expectedResult).withCoder(KvCoder.of(ReferenceShard.CODER, IterableCoder.of(new GATKReadCoder()))));
        DataflowTestUtils.keyIterableValueMatcher(grouped, pFinalExpected);

        p.run();
    }
}