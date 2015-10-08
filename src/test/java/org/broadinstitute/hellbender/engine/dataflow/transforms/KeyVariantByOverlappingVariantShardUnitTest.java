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
import org.broadinstitute.hellbender.engine.dataflow.coders.VariantCoder;
import org.broadinstitute.hellbender.engine.VariantShard;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

public final class KeyVariantByOverlappingVariantShardUnitTest extends BaseTest {

    @DataProvider(name = "keyedVariantShards")
    public Object[][] keyedVariantShards(){
        Object[][] data = new Object[2][];
        List<Class<?>> classes = Arrays.asList(Read.class, SAMRecord.class);
        for (int i = 0; i < classes.size(); ++i) {
            Class<?> c = classes.get(i);
            ReadsPreprocessingPipelineTestData testData = new ReadsPreprocessingPipelineTestData(c);

            List<Variant> variants = testData.getVariants();
            List<KV<VariantShard, Variant>> kvVariantShardVariant = testData.getKvVariantShardVariant();
            data[i] = new Object[]{variants, kvVariantShardVariant};
        }

        return data;
    }

    @Test(dataProvider = "keyedVariantShards")
    public void keyVariantsByVariantShardTest(List<Variant> variantList, List<KV<VariantShard, Variant>> kvVariantShardVariant) {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<Variant> pVariants = DataflowTestUtils.pCollectionCreateAndVerify(p, variantList, new VariantCoder());

        PCollection<KV<VariantShard, Variant>> kVariant = pVariants.apply(new KeyVariantByOverlappingVariantShard());
        DataflowAssert.that(kVariant).containsInAnyOrder(kvVariantShardVariant);
        p.run();
    }
}