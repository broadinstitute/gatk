package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.IterableCoder;
import com.google.cloud.dataflow.sdk.coders.KvCoder;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.collect.Maps;
import com.google.api.services.genomics.model.Read;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.engine.dataflow.ReadsPreprocessingPipelineTestData;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.engine.dataflow.coders.ReadContextDataCoder;
import org.broadinstitute.hellbender.engine.dataflow.coders.VariantCoder;
import org.broadinstitute.hellbender.engine.dataflow.datasources.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.FakeReferenceSource;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import static org.mockito.Mockito.*;

public final class AddContextDataToReadUnitTest extends BaseTest {

    @DataProvider(name = "bases")
    public Object[][] bases() {
        Object[][] data = new Object[2][];
        List<Class<?>> classes = Arrays.asList(Read.class, SAMRecord.class);
        for (int i = 0; i < classes.size(); ++i) {
            Class<?> c = classes.get(i);
            ReadsPreprocessingPipelineTestData testData = new ReadsPreprocessingPipelineTestData(c);

            List<GATKRead> reads = testData.getReads();
            List<KV<GATKRead, ReferenceBases>> kvReadRefBases = testData.getKvReadsRefBases();
            List<SimpleInterval> intervals = testData.getAllIntervals();
            List<Variant> variantList = testData.getVariants();
            List<KV<GATKRead, Iterable<Variant>>> kvReadiVariant = testData.getKvReadiVariant();
            List<KV<GATKRead, ReadContextData>> kvReadContextData = testData.getKvReadContextData();
            data[i] = new Object[]{reads, variantList, kvReadRefBases, kvReadContextData, intervals, kvReadiVariant};
        }
        return data;
    }

    @Test(dataProvider = "bases")
    public void addContextDataTest(List<GATKRead> reads, List<Variant> variantList,
                                   List<KV<GATKRead, ReferenceBases>> kvReadRefBases, List<KV<GATKRead, ReadContextData>> kvReadContextData,
                                   List<SimpleInterval> intervals, List<KV<GATKRead, Iterable<Variant>>> kvReadiVariant) {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<GATKRead> pReads = DataflowTestUtils.pCollectionCreateAndVerify(p, reads, new GATKReadCoder());
        PCollection<KV<GATKRead, ReferenceBases>> pReadRef = DataflowTestUtils.pCollectionCreateAndVerify(p, kvReadRefBases,
                KvCoder.of(new GATKReadCoder(), SerializableCoder.of(ReferenceBases.class)));

        PCollection<KV<GATKRead, Iterable<Variant>>> pReadVariants =
                p.apply(Create.of(kvReadiVariant).withCoder(KvCoder.of(new GATKReadCoder(), IterableCoder.of(new VariantCoder()))));

        PCollection<KV<GATKRead, ReadContextData>> joinedResults = AddContextDataToRead.join(pReads, pReadRef, pReadVariants);
        PCollection<KV<GATKRead, ReadContextData>> pkvReadContextData = p.apply(Create.of(kvReadContextData).withCoder(KvCoder.of(new GATKReadCoder(), new ReadContextDataCoder())));
        DataflowTestUtils.keyReadContextDataMatcher(joinedResults, pkvReadContextData);
        p.run();
    }

    @Test(dataProvider = "bases")
    public void fullTest(List<GATKRead> reads, List<Variant> variantList,
                      List<KV<GATKRead, ReferenceBases>> kvReadRefBases, List<KV<GATKRead, ReadContextData>> kvReadContextData,
                      List<SimpleInterval> intervals, List<KV<GATKRead, Iterable<Variant>>> kvReadiVariant) {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<GATKRead> pReads = DataflowTestUtils.pCollectionCreateAndVerify(p, reads, new GATKReadCoder());

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
        RefAPISource.setRefAPISource(mockSource);
        PCollection<KV<GATKRead, ReadContextData>> result = AddContextDataToRead.add(pReads, /*mockSource,*/ refAPIMetadata, mockVariantsSource);
        PCollection<KV<GATKRead, ReadContextData>> pkvReadContextData = p.apply(Create.of(kvReadContextData).withCoder(KvCoder.of(new GATKReadCoder(), new ReadContextDataCoder())));
        DataflowTestUtils.keyReadContextDataMatcher(result, pkvReadContextData);
        p.run();
    }
}