package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.collect.Maps;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestData;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.engine.dataflow.datasources.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;
import java.util.Map;

import static org.mockito.Matchers.any;
import static org.mockito.Matchers.eq;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;
import static org.mockito.Mockito.withSettings;

public final class RefBasesForReadsTest {
    @DataProvider(name = "bases")
    public Object[][] bases(){
        DataflowTestData testData = new DataflowTestData();

        List<GATKRead> reads = testData.getReads();
        List<SimpleInterval> intervals = testData.getAllIntervals();
        List<KV<GATKRead, ReferenceBases>> kvReadRefBases = testData.getKvReadsRefBases();

        return new Object[][]{
                {reads, kvReadRefBases, intervals},
        };
    }

    @Test(dataProvider = "bases")
    public void refBasesTest(List<GATKRead> reads, List<KV<GATKRead, ReferenceBases>> kvReadRefBases,
                             List<SimpleInterval> intervals) {
        Pipeline p = TestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<GATKRead> pReads = DataflowTestUtils.PCollectionCreateAndVerify(p, reads);

        RefAPISource mockSource = mock(RefAPISource.class, withSettings().serializable());
        for (SimpleInterval i : intervals) {
            when(mockSource.getReferenceBases(any(PipelineOptions.class), any(RefAPIMetadata.class), eq(i))).thenReturn(FakeReferenceSource.bases(i));
        }

        String referenceName = "refName";
        String crazyName = "0xbjfjd23f";
        Map<String, String> referenceNameToIdTable = Maps.newHashMap();
        referenceNameToIdTable.put(referenceName, crazyName);
        RefAPIMetadata refAPIMetadata = new RefAPIMetadata(referenceName, referenceNameToIdTable);
        PCollection<KV<GATKRead, ReferenceBases>> result = RefBasesForReads.addBases(pReads, mockSource, refAPIMetadata);
        DataflowAssert.that(result).containsInAnyOrder(kvReadRefBases);
        p.run();
    }

}