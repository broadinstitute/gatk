package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.IterableCoder;
import com.google.cloud.dataflow.sdk.coders.KvCoder;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.ReferenceShard;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.tools.ReadsPreprocessingPipelineTestData;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.FakeReferenceSource;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.mockito.Matchers.any;
import static org.mockito.Matchers.eq;
import static org.mockito.Mockito.*;

public final class RefBasesFromAPIUnitTest extends BaseTest {

    private ReferenceMultiSource createMockReferenceDataflowSource(final List<SimpleInterval> intervals,
                                                                      final SerializableFunction<GATKRead, SimpleInterval> referenceWindowingFunction) throws IOException {
        ReferenceMultiSource mockSource = mock(ReferenceMultiSource.class, withSettings().serializable());
        for (SimpleInterval interval : intervals) {
            when(mockSource.getReferenceBases(any(PipelineOptions.class), eq(interval))).thenReturn(FakeReferenceSource.bases(interval));
        }
        when(mockSource.getReferenceWindowFunction()).thenReturn(referenceWindowingFunction);
        return mockSource;
    }

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
                             List<SimpleInterval> intervals, List<KV<ReferenceBases, Iterable<GATKRead>>> kvRefBasesiReads) throws IOException {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        final ReferenceMultiSource mockSource = createMockReferenceDataflowSource(intervals, ReferenceWindowFunctions.IDENTITY_FUNCTION);

        PCollection<KV<ReferenceShard, Iterable<GATKRead>>> pInput = p.apply("pInput.Create", Create.of(kvRefShardiReads).withCoder(KvCoder.of(ReferenceShard.CODER, IterableCoder.of(new GATKReadCoder()))));

        PCollection<KV<ReferenceBases, Iterable<GATKRead>>> kvpCollection = RefBasesFromAPI.getBasesForShard(pInput, mockSource);
        PCollection<KV<ReferenceBases, Iterable<GATKRead>>> pkvRefBasesiReads = p.apply("pkvRefBasesiReads.Create", Create.of(kvRefBasesiReads).withCoder(KvCoder.of(SerializableCoder.of(ReferenceBases.class), IterableCoder.of(new GATKReadCoder()))));
        DataflowTestUtils.keyIterableValueMatcher(kvpCollection, pkvRefBasesiReads);

        p.run();
    }

    @Test(dataProvider = "bases")
    public void noReadsTest(List<KV<ReferenceShard, Iterable<GATKRead>>> kvRefShardiReads,
                             List<SimpleInterval> intervals, List<KV<ReferenceBases, Iterable<GATKRead>>> kvRefBasesiReads) throws IOException {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        final ReferenceMultiSource mockSource = createMockReferenceDataflowSource(intervals, ReferenceWindowFunctions.IDENTITY_FUNCTION);

        List<KV<ReferenceShard, Iterable<GATKRead>>> noReads = Arrays.asList(
                KV.of(new ReferenceShard(0, "1"), Lists.newArrayList()));

        PCollection<KV<ReferenceShard, Iterable<GATKRead>>> pInput = p.apply(Create.of(noReads).withCoder(KvCoder.of(ReferenceShard.CODER, IterableCoder.of(new GATKReadCoder()))));

        // We expect an exception to be thrown if there is a problem. If no error is thrown, it's fine.
        RefBasesFromAPI.getBasesForShard(pInput, mockSource);

        p.run();
    }

    @DataProvider(name = "RefBasesFromAPIWithCustomWindowFunctionTestData")
    public Object[][] refBasesFromAPIWithCustomWindowFunctionTestData() {
        final ReferenceShard shard = new ReferenceShard(0, "1");
        final List<GATKRead> samReads = Arrays.asList(
                ReadsPreprocessingPipelineTestData.makeRead("1", 500, 100, 0, SAMRecord.class),
                ReadsPreprocessingPipelineTestData.makeRead("1", 1000, 100, 1, SAMRecord.class),
                ReadsPreprocessingPipelineTestData.makeRead("1", 1050, 100, 2, SAMRecord.class)
        );
        final List<GATKRead> googleReads = Arrays.asList(
                ReadsPreprocessingPipelineTestData.makeRead("1", 500, 100, 0, Read.class),
                ReadsPreprocessingPipelineTestData.makeRead("1", 1000, 100, 1, Read.class),
                ReadsPreprocessingPipelineTestData.makeRead("1", 1050, 100, 2, Read.class)
        );

        final List<Object[]> testCases = new ArrayList<>();
        for ( final List<GATKRead> reads : Arrays.asList(samReads, googleReads) ) {
            // Test case layout: input shard + reads, reference window function to apply, expected reference interval post-transform

            // Identity function
            testCases.add( new Object[]{ KV.of(shard, reads), ReferenceWindowFunctions.IDENTITY_FUNCTION, new SimpleInterval("1", 500, 1149) });
            // Expand reads by 1 base on each side
            testCases.add( new Object[]{ KV.of(shard, reads), new ReferenceWindowFunctions.FixedWindowFunction(1, 1), new SimpleInterval("1", 499, 1150) });
            // Expand reads by 3 bases on the left and 5 bases on the right
            testCases.add( new Object[]{ KV.of(shard, reads), new ReferenceWindowFunctions.FixedWindowFunction(3, 5), new SimpleInterval("1", 497, 1154) });
            // Expand reads by 30 bases on the left and 10 bases on the right
            testCases.add( new Object[]{ KV.of(shard, reads), new ReferenceWindowFunctions.FixedWindowFunction(30, 10), new SimpleInterval("1", 470, 1159) });
            // Expand reads by 1000 bases on each side (goes off contig bounds)
            testCases.add( new Object[]{ KV.of(shard, reads), new ReferenceWindowFunctions.FixedWindowFunction(1000, 1000), new SimpleInterval("1", 1, 2149) });
        }
        return testCases.toArray(new Object[][]{});
    }

    @Test(dataProvider = "RefBasesFromAPIWithCustomWindowFunctionTestData")
    public void testRefBasesFromAPIWithCustomWindowFunction( final KV<ReferenceShard, Iterable<GATKRead>> inputShard, final SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction, final SimpleInterval expectedReferenceInterval ) throws IOException {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        // Assuming everything works properly, we only need a mock response to the expected reference interval,
        // since that is what should end up being queried in ReferenceDataflowSource. If we later try to query an
        // unexpected interval, we'll detect it via an Assert.
        final ReferenceMultiSource mockSource = createMockReferenceDataflowSource(Arrays.asList(expectedReferenceInterval), referenceWindowFunction);

        PCollection<KV<ReferenceShard, Iterable<GATKRead>>> pInput = p.apply("pInput.Create", Create.of(inputShard).withCoder(KvCoder.of(ReferenceShard.CODER, IterableCoder.of(new GATKReadCoder()))));
        PCollection<KV<ReferenceBases, Iterable<GATKRead>>> actualResult = RefBasesFromAPI.getBasesForShard(pInput, mockSource);

        DataflowAssert.that(actualResult).satisfies((Iterable<KV<ReferenceBases, Iterable<GATKRead>>> input) -> {
            for (KV<ReferenceBases, Iterable<GATKRead>> kvPair : input) {
                Assert.assertNotNull(kvPair.getKey(), "Null ReferenceBases in KV pair indicates that reference query in ReferenceDataflowSource used an unexpected/incorrect interval (mock ReferenceDataflowSource returned null)");
                Assert.assertEquals(kvPair.getKey().getInterval(), expectedReferenceInterval, "Wrong interval for ReferenceBases object after applying RefBasesFromAPI");
            }
            return null;
        });

        p.run();
    }
}
