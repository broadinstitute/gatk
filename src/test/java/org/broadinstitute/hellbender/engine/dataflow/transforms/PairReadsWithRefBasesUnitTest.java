package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.*;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.api.services.genomics.model.Read;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.engine.dataflow.ReadsPreprocessingPipelineTestData;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefWindowFunctions;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.test.BaseTest;

import static org.broadinstitute.hellbender.engine.dataflow.ReadsPreprocessingPipelineTestData.makeRead;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class PairReadsWithRefBasesUnitTest extends BaseTest {
    @DataProvider(name = "bases")
    public Object[][] bases(){
        Object[][] data = new Object[2][];
        List<Class<?>> classes = Arrays.asList(Read.class, SAMRecord.class);
        for (int i = 0; i < classes.size(); ++i) {
            Class<?> c = classes.get(i);
            ReadsPreprocessingPipelineTestData testData = new ReadsPreprocessingPipelineTestData(c);

            List<KV<ReferenceBases, Iterable<GATKRead>>> kvRefBasesiReads = testData.getKvRefBasesiReads();
            List<KV<GATKRead, ReferenceBases>> kvReadsRefBases = testData.getKvReadsRefBases();
            data[i] = new Object[]{kvRefBasesiReads, kvReadsRefBases};
        }
        return data;
    }

    @Test(dataProvider = "bases")
    public void refBasesTest(List<KV<ReferenceBases, Iterable<GATKRead>>> kvRefBasesiReads,
                             List<KV<GATKRead,ReferenceBases>> kvReadsRefBases) {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<KV<ReferenceBases, Iterable<GATKRead>>> pInput =
                p.apply(Create.of(kvRefBasesiReads).withCoder(
                        KvCoder.of(SerializableCoder.of(ReferenceBases.class), IterableCoder.of(new GATKReadCoder()))));

        PCollection<KV<GATKRead, ReferenceBases>> result = pInput.apply(new PairReadWithRefBases());
        DataflowAssert.that(result).containsInAnyOrder(kvReadsRefBases);

        p.run();
    }

    @DataProvider(name = "PairReadsWithRefBasesWithCustomWindowFunctionTestData")
    public Object[][] pairReadsWithRefBasesWithCustomWindowFunctionTestData() {
        // Using this made-up sequence as our reference for these tests.
        final byte[] bases = "AGCCTTTCGAACTGAGCCCGTTCCTGGGGTTATACCCGGCTTTTGGCGCT".getBytes();
        final ReferenceBases refBases = new ReferenceBases(bases, new SimpleInterval("1", 1, 50));

        final List<Object[]> testCases = new ArrayList<>();
        for ( final Class<?> readImplementation : Arrays.asList(SAMRecord.class, Read.class) ) {
            // Test case layout: Ref bases + read, reference window function to apply, expected ref bases for read, expected ref interval for read

            // Read at start of contig, identity function
            testCases.add(new Object[]{ KV.of(refBases, Arrays.asList(makeRead("1", 1, 10, 0, readImplementation))), RefWindowFunctions.IDENTITY_FUNCTION, "AGCCTTTCGA".getBytes(), new SimpleInterval("1", 1, 10) });
            // Read at start of contig, expand by 1 base on each side (goes off contig bounds)
            testCases.add(new Object[]{ KV.of(refBases, Arrays.asList(makeRead("1", 1, 10, 0, readImplementation))), new RefWindowFunctions.FixedWindowFunction(1, 1), "AGCCTTTCGAA".getBytes(), new SimpleInterval("1", 1, 11) });
            // Read at start of contig, expand by 3 bases on the left and 5 bases on the right (goes off contig bounds)
            testCases.add(new Object[]{ KV.of(refBases, Arrays.asList(makeRead("1", 1, 10, 0, readImplementation))), new RefWindowFunctions.FixedWindowFunction(3, 5), "AGCCTTTCGAACTGA".getBytes(), new SimpleInterval("1", 1, 15) });
            // Read in middle of contig, identity function
            testCases.add(new Object[]{ KV.of(refBases, Arrays.asList(makeRead("1", 20, 11, 0, readImplementation))), RefWindowFunctions.IDENTITY_FUNCTION, "GTTCCTGGGGT".getBytes(), new SimpleInterval("1", 20, 30) });
            // Read in middle of contig, expand by 1 base on each side
            testCases.add(new Object[]{ KV.of(refBases, Arrays.asList(makeRead("1", 20, 11, 0, readImplementation))), new RefWindowFunctions.FixedWindowFunction(1, 1), "CGTTCCTGGGGTT".getBytes(), new SimpleInterval("1", 19, 31) });
            // Read in middle of contig, expand by 3 bases on the left and 5 bases on the right
            testCases.add(new Object[]{ KV.of(refBases, Arrays.asList(makeRead("1", 20, 11, 0, readImplementation))), new RefWindowFunctions.FixedWindowFunction(3, 5), "CCCGTTCCTGGGGTTATAC".getBytes(), new SimpleInterval("1", 17, 35) });
            // Read in middle of contig, expand by 30 bases on the left and 10 bases on the right (goes off contig bounds)
            testCases.add(new Object[]{ KV.of(refBases, Arrays.asList(makeRead("1", 20, 11, 0, readImplementation))), new RefWindowFunctions.FixedWindowFunction(30, 10), "AGCCTTTCGAACTGAGCCCGTTCCTGGGGTTATACCCGGC".getBytes(), new SimpleInterval("1", 1, 40) });
        }
        return testCases.toArray(new Object[][]{});
    }

    @Test(dataProvider = "PairReadsWithRefBasesWithCustomWindowFunctionTestData")
    public void testPairReadsWithRefBasesWithCustomWindowFunction( final KV<ReferenceBases, Iterable<GATKRead>> input, final SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction, final byte[] expectedBases, final SimpleInterval expectedInterval ) {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<KV<ReferenceBases, Iterable<GATKRead>>> pInput = p.apply(Create.of(input)
                .withCoder(KvCoder.of(SerializableCoder.of(ReferenceBases.class), IterableCoder.of(new GATKReadCoder()))));

        PCollection<KV<GATKRead, ReferenceBases>> result = pInput.apply(new PairReadWithRefBases(referenceWindowFunction));

        DataflowAssert.that(result).satisfies((Iterable<KV<GATKRead, ReferenceBases>> resultElements) -> {
            for ( KV<GATKRead, ReferenceBases> kvPair : resultElements ) {
                Assert.assertNotNull(kvPair.getKey(), "Null read in transform result");
                Assert.assertNotNull(kvPair.getValue(), "Null ReferenceBases object paired with read");
                Assert.assertEquals(kvPair.getValue().getBases(), expectedBases, "Wrong bases in ReferenceBases object paired with read");
                Assert.assertEquals(kvPair.getValue().getInterval(), expectedInterval, "Wrong interval in ReferenceBases object paired with read");
            }
            return null;
        });

        p.run();
    }
}