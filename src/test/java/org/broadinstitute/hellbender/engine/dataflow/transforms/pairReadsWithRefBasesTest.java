package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.*;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestData;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;

public final class pairReadsWithRefBasesTest {
    @DataProvider(name = "bases")
    public Object[][] bases(){

        DataflowTestData testData = new DataflowTestData();

        List<KV<ReferenceBases, Iterable<GATKRead>>> kvRefBasesiReads = testData.getKvRefBasesiReads();
        List<KV<GATKRead, ReferenceBases>> kvReadsRefBases = testData.getKvReadsRefBases();

        return new Object[][]{
                {kvRefBasesiReads, kvReadsRefBases},
        };
    }

    @Test(dataProvider = "bases")
    public void refBasesTest(List<KV<ReferenceBases, Iterable<GATKRead>>> kvRefBasesiReads,
                             List<KV<GATKRead,ReferenceBases>> kvReadsRefBases) {
        Pipeline p = TestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<KV<ReferenceBases, Iterable<GATKRead>>> pInput =
                p.apply(Create.of(kvRefBasesiReads)).setCoder(
                        KvCoder.of(SerializableCoder.of(ReferenceBases.class), IterableCoder.of(new GATKReadCoder())));

        PCollection<KV<GATKRead, ReferenceBases>> result = pInput.apply(new PairReadWithRefBases());
        DataflowAssert.that(result).containsInAnyOrder(kvReadsRefBases);

        p.run();
    }

}