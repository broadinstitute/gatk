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
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;

public final class PairReadsAndVariantsTest {

    @DataProvider(name = "pairedReadsAndVariants")
    public Object[][] pairedReadsAndVariants(){
        DataflowTestData testData = new DataflowTestData();

        List<GATKRead> reads = testData.getReads();
        List<Variant> variantList = testData.getVariants();
        List<KV<GATKRead, Variant>> expected = testData.getKvReadVariant();

        return new Object[][]{
                {reads, variantList, expected},
        };
    }

    @Test(dataProvider = "pairedReadsAndVariants")
    public void pairReadsAndVariantsTest(List<GATKRead> reads, List<Variant> variantList, List<KV<GATKRead, Variant>> expected) {
        Pipeline p = TestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<GATKRead> pReads = DataflowTestUtils.PCollectionCreateAndVerify(p, reads);
        PCollection<Variant> pVariants = DataflowTestUtils.PCollectionCreateAndVerify(p, variantList);


        PCollection<KV<GATKRead, Variant>> readVariants = PairReadsAndVariants.pair(pReads, pVariants);
        DataflowAssert.that(readVariants).containsInAnyOrder(expected);

        p.run();
    }
}