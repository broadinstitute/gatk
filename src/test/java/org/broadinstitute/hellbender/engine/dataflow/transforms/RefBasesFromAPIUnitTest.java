package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.IterableCoder;
import com.google.cloud.dataflow.sdk.coders.KvCoder;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.engine.dataflow.ReadsPreprocessingPipelineTestData;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPIMetadata;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPISource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceShard;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.FakeReferenceSource;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import static org.mockito.Matchers.any;
import static org.mockito.Matchers.eq;
import static org.mockito.Mockito.*;

public final class RefBasesFromAPIUnitTest extends BaseTest {

    @DataProvider(name = "bases")
    public Object[][] bases(){
        Object[][] data = new Object[2][];
        List<Class<?>> classes = Arrays.asList(Read.class, SAMRecord.class);
        for (int i = 0; i < classes.size(); ++i) {
            Class<?> c = classes.get(i);
            ReadsPreprocessingPipelineTestData testData = new ReadsPreprocessingPipelineTestData(c);
            List<KV<ReferenceShard, Iterable<GATKRead>>> kvRefShardiReads = testData.getKvRefShardiReads();
            List<SimpleInterval> intervals = testData.getAllIntervals();
            List<KV<ReferenceBases, Iterable<GATKRead>>> kvRefBasesiReads = testData.getKvRefBasesiReads();
            data[i] = new Object[]{kvRefShardiReads, intervals, kvRefBasesiReads};
        }
        return data;
    }

    @Test(dataProvider = "bases")
    public void refBasesTest(List<KV<ReferenceShard, Iterable<GATKRead>>> kvRefShardiReads,
                             List<SimpleInterval> intervals, List<KV<ReferenceBases, Iterable<GATKRead>>> kvRefBasesiReads) {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        // Set up the mock for RefAPISource;
        RefAPISource mockSource = mock(RefAPISource.class, withSettings().serializable());
        for (SimpleInterval interval : intervals) {
            when(mockSource.getReferenceBases(any(PipelineOptions.class), any(RefAPIMetadata.class), eq(interval))).thenReturn(FakeReferenceSource.bases(interval));
        }

        String referenceName = "refName";
        String refId = "0xbjfjd23f";
        Map<String, String> referenceNameToIdTable = Maps.newHashMap();
        referenceNameToIdTable.put(referenceName, refId);
        RefAPIMetadata refAPIMetadata = new RefAPIMetadata(referenceName, referenceNameToIdTable);

        PCollection<KV<ReferenceShard, Iterable<GATKRead>>> pInput = p.apply("pInput.Create",Create.of(kvRefShardiReads).withCoder(KvCoder.of(ReferenceShard.CODER, IterableCoder.of(new GATKReadCoder()))));

        RefAPISource.setRefAPISource(mockSource);
        PCollection<KV<ReferenceBases, Iterable<GATKRead>>> kvpCollection = RefBasesFromAPI.getBasesForShard(pInput, refAPIMetadata);
        PCollection<KV<ReferenceBases, Iterable<GATKRead>>> pkvRefBasesiReads = p.apply("pkvRefBasesiReads.Create", Create.of(kvRefBasesiReads).withCoder(KvCoder.of(SerializableCoder.of(ReferenceBases.class), IterableCoder.of(new GATKReadCoder()))));
        DataflowTestUtils.keyIterableValueMatcher(kvpCollection, pkvRefBasesiReads);

        p.run();
    }
    @Test(dataProvider = "bases")
    public void noReadsTest(List<KV<ReferenceShard, Iterable<GATKRead>>> kvRefShardiReads,
                             List<SimpleInterval> intervals, List<KV<ReferenceBases, Iterable<GATKRead>>> kvRefBasesiReads) {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        // Set up the mock for RefAPISource;
        RefAPISource mockSource = mock(RefAPISource.class, withSettings().serializable());
        for (SimpleInterval interval : intervals) {
            when(mockSource.getReferenceBases(any(PipelineOptions.class), any(RefAPIMetadata.class), eq(interval))).thenReturn(FakeReferenceSource.bases(interval));
        }

        String referenceName = "refName";
        String refId = "0xbjfjd23f";
        Map<String, String> referenceNameToIdTable = Maps.newHashMap();
        referenceNameToIdTable.put(referenceName, refId);
        RefAPIMetadata refAPIMetadata = new RefAPIMetadata(referenceName, referenceNameToIdTable);


        List<KV<ReferenceShard, Iterable<GATKRead>>> noReads = Arrays.asList(
                KV.of(new ReferenceShard(0, "1"), Lists.newArrayList()));

        PCollection<KV<ReferenceShard, Iterable<GATKRead>>> pInput = p.apply(Create.of(noReads).withCoder(KvCoder.of(ReferenceShard.CODER, IterableCoder.of(new GATKReadCoder()))));

        RefAPISource.setRefAPISource(mockSource);
        // We expect an exception to be thrown if there is a problem. If no error is thrown, it's fine.
        RefBasesFromAPI.getBasesForShard(pInput, refAPIMetadata);

        p.run();
    }
}
