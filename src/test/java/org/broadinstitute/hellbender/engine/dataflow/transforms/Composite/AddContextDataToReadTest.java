package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.IterableCoder;
import com.google.cloud.dataflow.sdk.coders.KvCoder;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.collect.Maps;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestData;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.engine.dataflow.coders.VariantCoder;
import org.broadinstitute.hellbender.engine.dataflow.datasources.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;
import java.util.Map;

import static org.mockito.Mockito.*;

public final class AddContextDataToReadTest {

    @DataProvider(name = "bases")
    public Object[][] bases() {
        DataflowTestData testData = new DataflowTestData();

        List<GATKRead> reads = testData.getReads();
        List<KV<GATKRead, ReferenceBases>> kvReadRefBases = testData.getKvReadsRefBases();
        List<SimpleInterval> intervals = testData.getAllIntervals();
        List<Variant> variantList = testData.getVariants();
        List<KV<GATKRead, Iterable<Variant>>> kvReadiVariant = testData.getKvReadiVariant();
        List<KV<GATKRead, ReadContextData>> kvReadContextData = testData.getKvReadContextData();

        return new Object[][]{
                {reads, variantList, kvReadRefBases, kvReadContextData, intervals, kvReadiVariant},
        };
    }

    @Test(dataProvider = "bases")
    public void addContextDataTest(List<GATKRead> reads, List<Variant> variantList,
                                   List<KV<GATKRead, ReferenceBases>> kvReadRefBases, List<KV<GATKRead, ReadContextData>> kvReadContextData,
                                   List<SimpleInterval> intervals, List<KV<GATKRead, Iterable<Variant>>> kvReadiVariant) {
        Pipeline p = TestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<GATKRead> pReads = DataflowTestUtils.PCollectionCreateAndVerify(p, reads);
        PCollection<KV<GATKRead, ReferenceBases>> pReadRef = DataflowTestUtils.PCollectionCreateAndVerify(p, kvReadRefBases);

        PCollection<KV<GATKRead, Iterable<Variant>>> pReadVariants =
                p.apply(Create.of(kvReadiVariant)).setCoder(KvCoder.of(new GATKReadCoder(), IterableCoder.of(new VariantCoder())));

        PCollection<KV<GATKRead, ReadContextData>> joinedResults = AddContextDataToRead.Join(pReads, pReadRef, pReadVariants);
        DataflowAssert.that(joinedResults).containsInAnyOrder(kvReadContextData);
        p.run();
    }

    @Test(dataProvider = "bases")
    public void fullTest(List<GATKRead> reads, List<Variant> variantList,
                      List<KV<GATKRead, ReferenceBases>> kvReadRefBases, List<KV<GATKRead, ReadContextData>> kvReadContextData,
                      List<SimpleInterval> intervals, List<KV<GATKRead, Iterable<Variant>>> kvReadiVariant) {
        Pipeline p = TestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<GATKRead> pReads = DataflowTestUtils.PCollectionCreateAndVerify(p, reads);

        PCollection<Variant> pVariant = p.apply(Create.of(variantList));
        VariantsDataflowSource mockVariantsSource = mock(VariantsDataflowSource.class);

        when(mockVariantsSource.getAllVariants()).thenReturn(pVariant);

        RefAPISource mockSource = mock(RefAPISource.class, withSettings().serializable());
        for (SimpleInterval i : intervals) {
            when(mockSource.getReferenceBases(any(PipelineOptions.class), any(RefAPIMetadata.class), eq(i))).thenReturn(FakeReferenceSource.bases(i));
        }

        String referenceName = "refName";
        String refId = "0xbjfjd23f";
        Map<String, String> referenceNameToIdTable = Maps.newHashMap();
        referenceNameToIdTable.put(referenceName, refId);
        RefAPIMetadata refAPIMetadata = new RefAPIMetadata(referenceName, referenceNameToIdTable);
        PCollection<KV<GATKRead, ReadContextData>> result = AddContextDataToRead.Add(pReads, mockSource, refAPIMetadata, mockVariantsSource);
        DataflowAssert.that(result).containsInAnyOrder(kvReadContextData);
        p.run();
    }
}