package org.broadinstitute.hellbender.engine.spark;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMRecord;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.ReadContextData;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.FakeReferenceSource;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import static org.mockito.Matchers.any;
import static org.mockito.Matchers.eq;
import static org.mockito.Mockito.*;

public class AddContextDataToReadSparkUnitTest extends BaseTest {
    @DataProvider(name = "bases")
    public Object[][] bases() {
        List<Class<?>> classes = Arrays.asList(Read.class, SAMRecord.class);
        JoinStrategy[] strategies = JoinStrategy.values();
        Object[][] data = new Object[classes.size() * strategies.length][];
        for (int i = 0; i < classes.size(); ++i) {
            Class<?> c = classes.get(i);
            ReadsPreprocessingPipelineSparkTestData testData = new ReadsPreprocessingPipelineSparkTestData(c);

            List<GATKRead> reads = testData.getReads();
            List<SimpleInterval> intervals = testData.getAllIntervals();
            List<GATKVariant> variantList = testData.getVariants();
            List<KV<GATKRead, ReadContextData>> expectedReadContextData = testData.getKvReadContextData();
            for (int j = 0; j < strategies.length; j++) {
                data[i * strategies.length + j] = new Object[]{reads, variantList, expectedReadContextData, intervals, strategies[j]};
            }
        }
        return data;
    }

    @Test(dataProvider = "bases", groups = "spark")
    public void addContextDataTest(List<GATKRead> reads, List<GATKVariant> variantList,
                                   List<KV<GATKRead, ReadContextData>> expectedReadContextData,
                                   List<SimpleInterval> intervals, JoinStrategy joinStrategy) throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        JavaRDD<GATKRead> rddReads = ctx.parallelize(reads);
        JavaRDD<GATKVariant> rddVariants = ctx.parallelize(variantList);

        ReferenceMultiSource mockSource = mock(ReferenceMultiSource.class, withSettings().serializable());
        for (SimpleInterval i : intervals) {
            when(mockSource.getReferenceBases(any(PipelineOptions.class), eq(i))).thenReturn(FakeReferenceSource.bases(i));
        }
        when(mockSource.getReferenceWindowFunction()).thenReturn(ReferenceWindowFunctions.IDENTITY_FUNCTION);

        JavaPairRDD<GATKRead, ReadContextData> rddActual = AddContextDataToReadSpark.add(rddReads, mockSource, rddVariants, joinStrategy);
        Map<GATKRead, ReadContextData> actual = rddActual.collectAsMap();

        Assert.assertEquals(actual.size(), expectedReadContextData.size());
        for (KV<GATKRead, ReadContextData> kv : expectedReadContextData) {
            ReadContextData readContextData = actual.get(kv.getKey());
            Assert.assertNotNull(readContextData);
            Assert.assertTrue(CollectionUtils.isEqualCollection(Lists.newArrayList(readContextData.getOverlappingVariants()),
                    Lists.newArrayList(kv.getValue().getOverlappingVariants())));
            Assert.assertEquals(readContextData.getOverlappingReferenceBases(), kv.getValue().getOverlappingReferenceBases());
        }
    }
}