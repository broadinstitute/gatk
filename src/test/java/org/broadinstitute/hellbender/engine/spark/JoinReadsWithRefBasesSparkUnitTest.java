package org.broadinstitute.hellbender.engine.spark;

import htsjdk.samtools.SAMRecord;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.utils.KV;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.testutils.FakeReferenceSource;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import static org.mockito.Mockito.*;

public class JoinReadsWithRefBasesSparkUnitTest extends GATKBaseTest {
    @DataProvider(name = "bases")
    public Object[][] bases(){
        List<Class<?>> classes = Collections.singletonList(SAMRecord.class);
        Object[][] data = new Object[classes.size()][];
        for (int i = 0; i < classes.size(); ++i) {
            Class<?> c = classes.get(i);
            ReadsPreprocessingPipelineSparkTestData testData = new ReadsPreprocessingPipelineSparkTestData(c);

            List<GATKRead> reads = testData.getReads();
            List<SimpleInterval> intervals = testData.getAllIntervals();
            List<KV<GATKRead, ReferenceBases>> kvReadRefBases = testData.getKvReadsRefBases();
            data[i] = new Object[]{reads, kvReadRefBases, intervals};
        }
        return data;
    }

    @Test(dataProvider = "bases", groups = "spark")
    public void refBasesShuffleTest(List<GATKRead> reads, List<KV<GATKRead, ReferenceBases>> kvReadRefBases,
                                    List<SimpleInterval> intervals) throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        JavaRDD<GATKRead> rddReads = ctx.parallelize(reads);

        ReferenceMultiSparkSource mockSource = mock(ReferenceMultiSparkSource.class, withSettings().serializable());
        for (SimpleInterval i : intervals) {
            when(mockSource.getReferenceBases(eq(i))).thenReturn(FakeReferenceSource.bases(i));
        }
        when(mockSource.getReferenceWindowFunction()).thenReturn(ReferenceWindowFunctions.IDENTITY_FUNCTION);

        JavaPairRDD<GATKRead, ReferenceBases> rddResult = ShuffleJoinReadsWithRefBases.addBases(mockSource, rddReads);
        Map<GATKRead, ReferenceBases> result = rddResult.collectAsMap();

        for (KV<GATKRead, ReferenceBases> kv : kvReadRefBases) {
            ReferenceBases referenceBases = result.get(kv.getKey());
            Assert.assertNotNull(referenceBases);
            Assert.assertEquals(kv.getValue(),referenceBases);
        }
    }

    @Test(dataProvider = "bases", groups = "spark")
    public void refBasesBroadcastTest(List<GATKRead> reads, List<KV<GATKRead, ReferenceBases>> kvReadRefBases,
                                      List<SimpleInterval> intervals) throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        JavaRDD<GATKRead> rddReads = ctx.parallelize(reads);

        ReferenceMultiSparkSource mockSource = mock(ReferenceMultiSparkSource.class, withSettings().serializable());
        for (SimpleInterval i : intervals) {
            when(mockSource.getReferenceBases(eq(i))).thenReturn(FakeReferenceSource.bases(i));
        }
        when(mockSource.getReferenceWindowFunction()).thenReturn(ReferenceWindowFunctions.IDENTITY_FUNCTION);

        JavaPairRDD<GATKRead, ReferenceBases> rddResult = BroadcastJoinReadsWithRefBases.addBases(mockSource, rddReads);
        Map<GATKRead, ReferenceBases> result = rddResult.collectAsMap();

        for (KV<GATKRead, ReferenceBases> kv : kvReadRefBases) {
            ReferenceBases referenceBases = result.get(kv.getKey());
            Assert.assertNotNull(referenceBases);
            Assert.assertEquals(kv.getValue(),referenceBases);
        }
    }
}