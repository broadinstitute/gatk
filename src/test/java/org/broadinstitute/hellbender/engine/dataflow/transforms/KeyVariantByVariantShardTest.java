package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestData;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.engine.dataflow.datasources.VariantShard;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;

public final class KeyVariantByVariantShardTest {

    @DataProvider(name = "keyedVariantShards")
    public Object[][] keyedVariantShards(){
        DataflowTestData testData = new DataflowTestData();

        List<Variant> variants = testData.getVariants();
        List<KV<VariantShard, Variant>> kvVariantShardVariant = testData.getKvVariantShardVariant();

        return new Object[][]{
                {variants, kvVariantShardVariant},
        };
    }

    @Test(dataProvider = "keyedVariantShards")
    public void keyVariantsByVariantShardTest(List<Variant> variantList, List<KV<VariantShard, Variant>> kvVariantShardVariant) {
        Pipeline p = TestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<Variant> pVariants = DataflowTestUtils.PCollectionCreateAndVerify(p, variantList);

        PCollection<KV<VariantShard, Variant>> kVariant = pVariants.apply(new KeyVariantByVariantShard());
        DataflowAssert.that(kVariant).containsInAnyOrder(kvVariantShardVariant);
        p.run();
    }
}