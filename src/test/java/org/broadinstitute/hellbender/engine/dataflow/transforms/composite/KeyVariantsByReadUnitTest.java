package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

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
import org.broadinstitute.hellbender.tools.ReadsPreprocessingPipelineTestData;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.engine.dataflow.coders.VariantCoder;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

public final class KeyVariantsByReadUnitTest extends BaseTest {

    @DataProvider(name = "variantsAndReads")
    public Object[][] variantsAndReads(){
        Object[][] data = new Object[2][];
        List<Class<?>> classes = Arrays.asList(Read.class, SAMRecord.class);
        for (int i = 0; i < classes.size(); ++i) {
            Class<?> c = classes.get(i);
            ReadsPreprocessingPipelineTestData testData = new ReadsPreprocessingPipelineTestData(c);

            List<GATKRead> reads = testData.getReads();
            List<Variant> variantList = testData.getVariants();
            List<KV<GATKRead, Iterable<Variant>>> kvReadiVariant = testData.getKvReadiVariantBroken();

            data[i] = new Object[]{reads, variantList, kvReadiVariant};
        }
        return data;
    }

    @Test(dataProvider = "variantsAndReads")
    public void keyVariantsByReadTest(List<GATKRead> reads, List<Variant> variantList,
                                   List<KV<GATKRead, Iterable<Variant>>> kvReadiVariant) {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<GATKRead> pReads = DataflowTestUtils.pCollectionCreateAndVerify(p, reads, new GATKReadCoder());
        PCollection<Variant> pVariant = DataflowTestUtils.pCollectionCreateAndVerify(p, variantList, new VariantCoder());

        PCollection<KV<GATKRead, Iterable<Variant>>> result = KeyVariantsByRead.key(pVariant, pReads);
        PCollection<KV<GATKRead, Iterable<Variant>>> pFinalExpected = p.apply(Create.of(kvReadiVariant).withCoder(KvCoder.of(new GATKReadCoder(), IterableCoder.of(new VariantCoder()))));
        DataflowTestUtils.keyIterableValueMatcher(result, pFinalExpected);

        p.run();
    }
}