package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.IterableCoder;
import com.google.cloud.dataflow.sdk.coders.KvCoder;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.engine.dataflow.ReadsPreprocessingPipelineTestData;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefWindowFunctions;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceShard;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class KeyReadsByRefShardUnitTest extends BaseTest {

    @DataProvider(name = "refShards")
    public Object[][] refShards(){
        Object[][] data = new Object[2][];
        List<Class<?>> classes = Arrays.asList(Read.class, SAMRecord.class);
        for (int i = 0; i < classes.size(); ++i) {
            Class<?> c = classes.get(i);
            ReadsPreprocessingPipelineTestData testData = new ReadsPreprocessingPipelineTestData(c);

            List<GATKRead> inputs = testData.getReads();
            List<KV<ReferenceShard, Iterable<GATKRead>>> kvs = testData.getKvRefShardiReads();
            data[i] = new Object[]{inputs, kvs};

        }
        return data;
    }

    @Test(dataProvider = "refShards")
    public void groupReadsForRefTest(List<GATKRead> reads, List<KV<ReferenceShard, Iterable<GATKRead>>> expectedResult) {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<GATKRead> pReads = DataflowTestUtils.pCollectionCreateAndVerify(p, reads, new GATKReadCoder());
        PCollection<KV<ReferenceShard, Iterable<GATKRead>>> grouped = pReads.apply(new KeyReadsByRefShard());

        PCollection<KV<ReferenceShard, Iterable<GATKRead>>> pFinalExpected = p.apply(Create.of(expectedResult).withCoder(KvCoder.of(ReferenceShard.CODER, IterableCoder.of(new GATKReadCoder()))));
        DataflowTestUtils.keyIterableValueMatcher(grouped, pFinalExpected);

        p.run();
    }

    @DataProvider(name = "KeyReadsByRefShardWithCustomWindowFunctionTestData")
    public Object[][] keyReadsByRefShardWithCustomWindowFunctionTestData() {
        final List<GATKRead> samReads = ReadsPreprocessingPipelineTestData.makeReferenceShardBoundaryReads(2, 3, SAMRecord.class);
        final List<GATKRead> googleReads = ReadsPreprocessingPipelineTestData.makeReferenceShardBoundaryReads(2, 3, Read.class);

        return new Object[][] {
                // Test case layout: reads, reference window function to apply, expected shard + reads pairs

                // Identity function, SAM reads
                { samReads, RefWindowFunctions.IDENTITY_FUNCTION, generateExpectedCustomWindowResult(samReads, 0, 0) },
                // Identity function, google reads
                { googleReads, RefWindowFunctions.IDENTITY_FUNCTION, generateExpectedCustomWindowResult(googleReads, 0, 0) },
                // Expand reads by 1 base on each side, SAM reads
                { samReads, new RefWindowFunctions.FixedWindowFunction(1, 1), generateExpectedCustomWindowResult(samReads, 1, 1) },
                // Expand reads by 1 base on each side, google reads
                { googleReads, new RefWindowFunctions.FixedWindowFunction(1, 1), generateExpectedCustomWindowResult(googleReads, 1, 1) },
                // Expand reads by 3 bases on the left and 5 bases on the right, SAM reads
                { samReads, new RefWindowFunctions.FixedWindowFunction(3, 5), generateExpectedCustomWindowResult(samReads, 3, 5) },
                // Expand reads by 3 bases on the left and 5 bases on the right, google reads
                { googleReads, new RefWindowFunctions.FixedWindowFunction(3, 5), generateExpectedCustomWindowResult(googleReads, 3, 5) },
        };
    }

    private List<KV<ReferenceShard, Iterable<GATKRead>>> generateExpectedCustomWindowResult( final List<GATKRead> reads, final int windowLeadingBases, final int windowTrailingBases ) {
        final Map<ReferenceShard, List<GATKRead>> shardMap = new HashMap<>();

        for ( final GATKRead read : reads ) {
            final SimpleInterval expandedReadInterval = new SimpleInterval(read.getContig(), Math.max(read.getStart() - windowLeadingBases, 1), read.getEnd() + windowTrailingBases);
            final ReferenceShard readShard = ReferenceShard.getShardNumberFromInterval(expandedReadInterval);

            if ( shardMap.containsKey(readShard) ) {
                shardMap.get(readShard).add(read);
            }
            else {
                shardMap.put(readShard, new ArrayList<>(Arrays.asList(read)));
            }
        }

        List<KV<ReferenceShard, Iterable<GATKRead>>> expectedResult = new ArrayList<>();
        for ( Map.Entry<ReferenceShard, List<GATKRead>> entry : shardMap.entrySet() ) {
            expectedResult.add(KV.of(entry.getKey(), entry.getValue()));
        }
        return expectedResult;
    }

    @Test(dataProvider = "KeyReadsByRefShardWithCustomWindowFunctionTestData")
    public void testKeyReadsByRefShardWithCustomWindowFunction(final List<GATKRead> reads, final SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction, final List<KV<ReferenceShard, Iterable<GATKRead>>> expectedResult ) {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<GATKRead> pReads = DataflowTestUtils.pCollectionCreateAndVerify(p, reads, new GATKReadCoder());
        PCollection<KV<ReferenceShard, Iterable<GATKRead>>> shardedReads = pReads.apply(new KeyReadsByRefShard(referenceWindowFunction));

        PCollection<KV<ReferenceShard, Iterable<GATKRead>>> pExpected = p.apply(Create.of(expectedResult).withCoder(KvCoder.of(ReferenceShard.CODER, IterableCoder.of(new GATKReadCoder()))));
        DataflowTestUtils.keyIterableValueMatcher(shardedReads, pExpected);

        p.run();
    }
}