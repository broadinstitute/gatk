package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

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
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

import static org.broadinstitute.hellbender.engine.dataflow.ReadsPreprocessingPipelineTestData.makeRead;
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

        final RefAPIMetadata refAPIMetadata = ReadsPreprocessingPipelineTestData.createRefAPIMetadata();
        RefAPISource.setRefAPISource(mockSource);
        PCollection<KV<GATKRead, ReadContextData>> result = AddContextDataToRead.add(pReads, /*mockSource,*/ new ReferenceDataflowSource(refAPIMetadata), mockVariantsSource);
        PCollection<KV<GATKRead, ReadContextData>> pkvReadContextData = p.apply(Create.of(kvReadContextData).withCoder(KvCoder.of(new GATKReadCoder(), new ReadContextDataCoder())));
        DataflowTestUtils.keyReadContextDataMatcher(result, pkvReadContextData);
        p.run();
    }

    @DataProvider(name = "AddContextDataWithCustomReferenceWindowFunctionTestData")
    public Object[][] addContextDataWithCustomReferenceWindowFunctionTestData() {
        // Using this made-up sequence as our reference for these tests.
        final String referenceBasesString = "AGCCTTTCGAACTGAGCCCGTTCCTGGGGTTATACCCGGCTTTTGGCGCT";
        final RefAPISource mockReferenceSource = mock(RefAPISource.class, withSettings().serializable());

        // Set up mock reference source to handle expected queries in our test cases below
        // (will return null for unexpected queries)
        for ( final SimpleInterval queryInterval : Arrays.asList(new SimpleInterval("1", 1, 10), new SimpleInterval("1", 1, 11), new SimpleInterval("1", 1, 15), new SimpleInterval("1", 20, 30), new SimpleInterval("1", 19, 31), new SimpleInterval("1", 17, 35), new SimpleInterval("1", 1, 40)) ) {
            when(mockReferenceSource.getReferenceBases(any(PipelineOptions.class), any(RefAPIMetadata.class), eq(queryInterval))).thenReturn(new ReferenceBases(referenceBasesString.substring(queryInterval.getStart() - 1, queryInterval.getEnd()).getBytes(), queryInterval));
        }

        final List<Object[]> testCases = new ArrayList<>();
        for ( final Class<?> readImplementation : Arrays.asList(SAMRecord.class, Read.class) ) {
            // Test case layout: read, mock reference source, reference window function to apply, expected ReferenceBases for read

            // Read at start of contig, identity function
            testCases.add(new Object[]{ makeRead("1", 1, 10, 0, readImplementation), mockReferenceSource, RefWindowFunctions.IDENTITY_FUNCTION, new ReferenceBases("AGCCTTTCGA".getBytes(), new SimpleInterval("1", 1, 10)) });
            // Read at start of contig, expand by 1 base on each side (goes off contig bounds)
            testCases.add(new Object[]{ makeRead("1", 1, 10, 0, readImplementation), mockReferenceSource, new RefWindowFunctions.FixedWindowFunction(1, 1), new ReferenceBases("AGCCTTTCGAA".getBytes(), new SimpleInterval("1", 1, 11)) });
            // Read at start of contig, expand by 3 bases on the left and 5 bases on the right (goes off contig bounds)
            testCases.add(new Object[]{ makeRead("1", 1, 10, 0, readImplementation), mockReferenceSource, new RefWindowFunctions.FixedWindowFunction(3, 5), new ReferenceBases("AGCCTTTCGAACTGA".getBytes(), new SimpleInterval("1", 1, 15)) });
            // Read in middle of contig, identity function
            testCases.add(new Object[]{ makeRead("1", 20, 11, 0, readImplementation), mockReferenceSource, RefWindowFunctions.IDENTITY_FUNCTION, new ReferenceBases("GTTCCTGGGGT".getBytes(), new SimpleInterval("1", 20, 30)) });
            // Read in middle of contig, expand by 1 base on each side
            testCases.add(new Object[]{ makeRead("1", 20, 11, 0, readImplementation), mockReferenceSource, new RefWindowFunctions.FixedWindowFunction(1, 1), new ReferenceBases("CGTTCCTGGGGTT".getBytes(), new SimpleInterval("1", 19, 31)) });
            // Read in middle of contig, expand by 3 bases on the left and 5 bases on the right
            testCases.add(new Object[]{ makeRead("1", 20, 11, 0, readImplementation), mockReferenceSource, new RefWindowFunctions.FixedWindowFunction(3, 5), new ReferenceBases("CCCGTTCCTGGGGTTATAC".getBytes(), new SimpleInterval("1", 17, 35)) });
            // Read in middle of contig, expand by 30 bases on the left and 10 bases on the right (goes off contig bounds)
            testCases.add(new Object[]{ makeRead("1", 20, 11, 0, readImplementation), mockReferenceSource, new RefWindowFunctions.FixedWindowFunction(30, 10), new ReferenceBases("AGCCTTTCGAACTGAGCCCGTTCCTGGGGTTATACCCGGC".getBytes(), new SimpleInterval("1", 1, 40)) });
        }
        return testCases.toArray(new Object[][]{});
    }

    @Test(dataProvider = "AddContextDataWithCustomReferenceWindowFunctionTestData")
    public void testAddContextDataWithCustomReferenceWindowFunction( final GATKRead read, final RefAPISource mockRefSource, final SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction, final ReferenceBases expectedReferenceBases ) {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<GATKRead> pReads = DataflowTestUtils.pCollectionCreateAndVerify(p, Arrays.asList(read), new GATKReadCoder());

        PCollection<Variant> pVariant = p.apply(Create.of(Collections.<Variant>emptyList()));
        VariantsDataflowSource mockVariantsSource = mock(VariantsDataflowSource.class);
        when(mockVariantsSource.getAllVariants()).thenReturn(pVariant);

        final RefAPIMetadata refAPIMetadata = ReadsPreprocessingPipelineTestData.createRefAPIMetadata(referenceWindowFunction);
        RefAPISource.setRefAPISource(mockRefSource);

        PCollection<KV<GATKRead, ReadContextData>> result = AddContextDataToRead.add(pReads, new ReferenceDataflowSource(refAPIMetadata), mockVariantsSource);

        DataflowAssert.that(result).satisfies((Iterable<KV<GATKRead, ReadContextData>> resultElements) -> {
            for ( KV<GATKRead, ReadContextData> kvPair : resultElements ) {
                Assert.assertNotNull(kvPair.getKey(), "Null read in KV pair after AddContextDataToRead");
                Assert.assertNotNull(kvPair.getValue(), "Null ReadContextData paired with read after AddContextDataToRead");
                Assert.assertNotNull(kvPair.getValue().getOverlappingReferenceBases(), "Null ReferenceBases in KV pair indicates that reference query in RefBaseFromAPI used an unexpected/incorrect interval (mock RefAPISource returned null)");
                Assert.assertEquals(kvPair.getValue().getOverlappingReferenceBases(), expectedReferenceBases, "Wrong ReferenceBases paired with read after AddContextDataToRead");
            }
            return null;
        });

        p.run();

    }
}