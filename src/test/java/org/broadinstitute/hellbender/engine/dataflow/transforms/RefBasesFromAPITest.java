package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.CollectionCoder;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestData;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.engine.dataflow.datasources.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import static org.mockito.Matchers.any;
import static org.mockito.Matchers.eq;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;
import static org.mockito.Mockito.withSettings;

public final class RefBasesFromAPITest {

    @DataProvider(name = "bases")
    public Object[][] bases(){
        DataflowTestData testData = new DataflowTestData();
        List<KV<ReferenceShard, Iterable<GATKRead>>> kvRefShardiReads = testData.getKvRefShardiReads();
        List<SimpleInterval> intervals = testData.getAllIntervals();
        List<KV<ReferenceBases, Iterable<GATKRead>>> kvRefBasesiReads = testData.getKvRefBasesiReads();

        return new Object[][]{
                {kvRefShardiReads, intervals, kvRefBasesiReads},
        };
    }

    @Test(dataProvider = "bases")
    public void refBasesTest(List<KV<ReferenceShard, Iterable<GATKRead>>> kvRefShardiReads,
                             List<SimpleInterval> intervals, List<KV<ReferenceBases, Iterable<GATKRead>>> kvRefBasesiReads) {
        Pipeline p = TestPipeline.create();
        p.getCoderRegistry().registerCoder(ArrayList.class, CollectionCoder.of(new GATKReadCoder()));
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

        PCollection<KV<ReferenceShard, Iterable<GATKRead>>> pInput = p.apply(Create.of(kvRefShardiReads));

        PCollection<KV<ReferenceBases, Iterable<GATKRead>>> kvpCollection = RefBasesFromAPI.getBasesForShard(pInput, mockSource, refAPIMetadata);
        DataflowAssert.that(kvpCollection).containsInAnyOrder(kvRefBasesiReads);

        p.run();
    }
    @Test(dataProvider = "bases")
    public void noReadsTest(List<KV<ReferenceShard, Iterable<GATKRead>>> kvRefShardiReads,
                             List<SimpleInterval> intervals, List<KV<ReferenceBases, Iterable<GATKRead>>> kvRefBasesiReads) {
        Pipeline p = TestPipeline.create();
        p.getCoderRegistry().registerCoder(ArrayList.class, CollectionCoder.of(new GATKReadCoder()));
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

        PCollection<KV<ReferenceShard, Iterable<GATKRead>>> pInput = p.apply(Create.of(noReads));

        // We expect an exception to be thrown if there is a problem. If no error is thrown, it's fine.
        RefBasesFromAPI.getBasesForShard(pInput, mockSource, refAPIMetadata);

        p.run();
    }
}
