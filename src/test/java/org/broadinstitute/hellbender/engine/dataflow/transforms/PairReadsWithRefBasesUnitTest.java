package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.*;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.api.services.genomics.model.Read;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.engine.dataflow.ReadsPreprocessingPipelineTestData;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

public final class PairReadsWithRefBasesUnitTest extends BaseTest {
    @DataProvider(name = "bases")
    public Object[][] bases(){
        Object[][] data = new Object[2][];
        List<Class<?>> classes = Arrays.asList(Read.class, SAMRecord.class);
        for (int i = 0; i < classes.size(); ++i) {
            Class<?> c = classes.get(i);
            ReadsPreprocessingPipelineTestData testData = new ReadsPreprocessingPipelineTestData(c);

            List<KV<ReferenceBases, Iterable<GATKRead>>> kvRefBasesiReads = testData.getKvRefBasesiReads();
            List<KV<GATKRead, ReferenceBases>> kvReadsRefBases = testData.getKvReadsRefBases();
            data[i] = new Object[]{kvRefBasesiReads, kvReadsRefBases};
        }
        return data;
    }

    @Test(dataProvider = "bases")
    public void refBasesTest(List<KV<ReferenceBases, Iterable<GATKRead>>> kvRefBasesiReads,
                             List<KV<GATKRead,ReferenceBases>> kvReadsRefBases) {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<KV<ReferenceBases, Iterable<GATKRead>>> pInput =
                p.apply(Create.of(kvRefBasesiReads).withCoder(
                        KvCoder.of(SerializableCoder.of(ReferenceBases.class), IterableCoder.of(new GATKReadCoder()))));

        PCollection<KV<GATKRead, ReferenceBases>> result = pInput.apply(new PairReadWithRefBases());
        DataflowAssert.that(result).containsInAnyOrder(kvReadsRefBases);

        p.run();
    }

}