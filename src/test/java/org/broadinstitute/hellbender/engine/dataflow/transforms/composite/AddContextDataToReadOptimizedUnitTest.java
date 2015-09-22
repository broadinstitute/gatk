package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

import static org.mockito.Matchers.any;
import static org.mockito.Matchers.eq;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;
import static org.mockito.Mockito.withSettings;
import static org.testng.Assert.*;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.transforms.Count;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.common.collect.Lists;

import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefWindowFunctions;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceDataflowSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.FakeReferenceSource;
import org.broadinstitute.hellbender.utils.variant.SkeletonVariant;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.UUID;

public class AddContextDataToReadOptimizedUnitTest extends BaseTest implements Serializable{
    private static final long serialVersionUID = 1L;


    // "cloud" because subdivideAndFillReads gets the reads from the cloud.
    @Test(groups = {"cloud"})
    public void testSubdivideAndFillReads() throws IOException {
        String bam = "gs://hellbender/test/resources/org/broadinstitute/hellbender/tools/BQSR/NA12878.chr17_69k_70k.dictFix.bam";
        int outputShardSize = 50;
        int margin = 1000;
        AddContextDataToReadOptimized.ContextShard shard = new AddContextDataToReadOptimized.ContextShard(new SimpleInterval("17",69000,69100));
        List<AddContextDataToReadOptimized.ContextShard> shards = new ArrayList<>();
        shards.add(shard);

        // 809R9ABXX101220:5:46:2178:137187 at 69032
        // 809R9ABXX101220:5:61:7847:167387 at 69038
        // then nothing until 69107
        HashSet<String> expectedReads = new HashSet<>();
        expectedReads.add("809R9ABXX101220:5:46:2178:137187");
        expectedReads.add("809R9ABXX101220:5:61:7847:167387");


        Pipeline p = GATKTestPipeline.create();
        p.getOptions().as(GCSOptions.class).setApiKey(System.getenv("HELLBENDER_TEST_APIKEY"));
        DataflowUtils.registerGATKCoders(p);
        PCollection<AddContextDataToReadOptimized.ContextShard> pShards = p.apply(Create.of(shards));
        pShards .apply(ParDo.of(AddContextDataToReadOptimized.subdivideAndFillReads(bam, outputShardSize, margin, null)))
                .apply(ParDo.of(new DoFn<AddContextDataToReadOptimized.ContextShard, String>() {
                    private static final long serialVersionUID = 1L;
                    @Override
                    public void processElement(ProcessContext c) throws Exception {
                        for (GATKRead r : c.element().reads) {
                            String n = r.getName();
                            System.out.println("Found "+n);
                            assertTrue(expectedReads.contains(n));
                            expectedReads.remove(n);
                            c.output(n);
                        }
                    }
                }))
                // we found all the reads we wanted
                .apply(Count.globally())
                .apply(ParDo.of(new DoFn<Long, Void>() {
                    private static final long serialVersionUID = 1L;
                                    @Override
                                    public void processElement(ProcessContext c) throws Exception {
                                        assertEquals(c.element().longValue(), 2L);
                                    }
                                }));
        p.run();
    }

    @Test
    public void testFillContext() throws IOException {
        GATKRead read = ArtificialReadUtils.createSamBackedReadWithUUID(new UUID(1, 2), "MADE:UP:READ", "17", 69000, 3);
        ArrayList<Variant> variants = makeTestVariant(69001, 69001);
        AddContextDataToReadOptimized.ContextShard shard = new AddContextDataToReadOptimized.ContextShard(new SimpleInterval("17",69000,69100))
                .withVariants(variants)
        .withReads(Lists.newArrayList(read));
        List < AddContextDataToReadOptimized.ContextShard> shards = new ArrayList<>();
        shards.add(shard);

        Pipeline p = GATKTestPipeline.create();
        p.getOptions().as(GCSOptions.class).setApiKey(System.getenv("HELLBENDER_TEST_APIKEY"));
        DataflowUtils.registerGATKCoders(p);
        PCollection<AddContextDataToReadOptimized.ContextShard> pShards = p.apply(Create.of(shards));

        ReferenceDataflowSource mockSource = mock(ReferenceDataflowSource.class, withSettings().serializable());
        SimpleInterval refInterval =  new SimpleInterval("17",69000,69002);
        when(mockSource.getReferenceBases(any(PipelineOptions.class), eq(refInterval))).thenReturn(FakeReferenceSource.bases(refInterval));
        when(mockSource.getReferenceWindowFunction()).thenReturn(RefWindowFunctions.IDENTITY_FUNCTION);

        pShards.apply(ParDo.of(AddContextDataToReadOptimized.fillContext(mockSource)))
                .apply(ParDo.of(new DoFn<AddContextDataToReadOptimized.ContextShard, Void>() {
                    private static final long serialVersionUID = 1L;
                    @Override
                    public void processElement(ProcessContext c) throws Exception {
                        ArrayList<ReadContextData> readContext = c.element().readContext;
                        assertEquals(readContext.get(0).getOverlappingReferenceBases().getInterval(), refInterval);
                        assertEquals(readContext.get(0).getOverlappingVariants(), variants);
                    }
                }));

        p.run();
    }


    @DataProvider(name="fillVariants")
    public Object[][] fillVariantsData() {
        List<SimpleInterval> intervals = Lists.newArrayList(
                new SimpleInterval("17", 100, 200)
        );
        ArrayList<Variant> variants1 = makeTestVariant(170, 180);
        List<AddContextDataToReadOptimized.ContextShard> expectedShards = new ArrayList<>();
        AddContextDataToReadOptimized.ContextShard shard1 = new AddContextDataToReadOptimized.ContextShard(intervals.get(0)).withVariants(variants1);
        expectedShards.add(shard1);

        ArrayList<Variant> variants2 = makeTestVariant(270, 280);
        List<AddContextDataToReadOptimized.ContextShard> expectedShards2 = new ArrayList<>();
        AddContextDataToReadOptimized.ContextShard shard2 = new AddContextDataToReadOptimized.ContextShard(intervals.get(0)).withVariants(variants2);
        expectedShards2.add(shard2);

        List<AddContextDataToReadOptimized.ContextShard> emptyShard = new ArrayList<>();
        AddContextDataToReadOptimized.ContextShard shard3 = new AddContextDataToReadOptimized.ContextShard(intervals.get(0)).withVariants(new ArrayList<>());
        emptyShard.add(shard3);

        return new Object[][] {
                { intervals, variants1, 100, expectedShards },
                { intervals, variants2, 100, expectedShards2 },
                { intervals, variants2,  50, emptyShard }
        };
    }

    private ArrayList<Variant> makeTestVariant(int begin, int end) {
        return Lists.newArrayList(
                new SkeletonVariant(new SimpleInterval("17", begin, end), true, false, new UUID(begin, end))
        );
    }

    @Test(dataProvider="fillVariants")
    public void testFillVariants(List<SimpleInterval> intervals, ArrayList<Variant> variants, int margin, List<AddContextDataToReadOptimized.ContextShard> expectedShards) {
        ArrayList<AddContextDataToReadOptimized.ContextShard> shards = AddContextDataToReadOptimized.fillVariants(intervals, variants, margin);
        assertEquals(shards.size(), expectedShards.size());
        for (int i=0; i<shards.size(); i++) {
            AddContextDataToReadOptimized.ContextShard a = shards.get(i);
            AddContextDataToReadOptimized.ContextShard b = expectedShards.get(i);
            assertEquals(a.interval, b.interval);
            assertEquals(a.variants.getOverlapping(a.interval), b.variants.getOverlapping(a.interval));
            assertEquals(a.reads, b.reads);
            assertEquals(a.readContext, b.readContext);
        }
    }

}