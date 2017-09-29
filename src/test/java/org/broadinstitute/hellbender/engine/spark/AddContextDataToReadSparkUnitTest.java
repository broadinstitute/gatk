package org.broadinstitute.hellbender.engine.spark;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.ReadContextData;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.FakeReferenceSource;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.mockito.invocation.InvocationOnMock;
import org.mockito.stubbing.Answer;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.io.Serializable;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import static org.mockito.Matchers.any;
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
            List<GATKVariant> variantList = testData.getVariants();
            List<KV<GATKRead, ReadContextData>> expectedReadContextData = testData.getKvReadContextData();
            for (int j = 0; j < strategies.length; j++) {
                data[i * strategies.length + j] = new Object[]{reads, variantList, expectedReadContextData, strategies[j]};
            }
        }
        return data;
    }

    @Test(dataProvider = "bases", groups = "spark")
    public void addContextDataTest(List<GATKRead> reads, List<GATKVariant> variantList,
                                   List<KV<GATKRead, ReadContextData>> expectedReadContextData,
                                   JoinStrategy joinStrategy) throws IOException {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        JavaRDD<GATKRead> rddReads = ctx.parallelize(reads);
        JavaRDD<GATKVariant> rddVariants = ctx.parallelize(variantList);

        ReferenceMultiSource mockSource = mock(ReferenceMultiSource.class, withSettings().serializable());
        when(mockSource.getReferenceBases(any(PipelineOptions.class), any())).then(new ReferenceBasesAnswer());
        when(mockSource.getReferenceWindowFunction()).thenReturn(ReferenceWindowFunctions.IDENTITY_FUNCTION);
        SAMSequenceDictionary sd = new SAMSequenceDictionary(Lists.newArrayList(new SAMSequenceRecord("1", 100000), new SAMSequenceRecord("2", 100000)));
        when(mockSource.getReferenceSequenceDictionary(null)).thenReturn(sd);

        JavaPairRDD<GATKRead, ReadContextData> rddActual = AddContextDataToReadSpark.add(ctx, rddReads, mockSource, rddVariants, null, joinStrategy,
                sd, 10000, 1000);
        Map<GATKRead, ReadContextData> actual = rddActual.collectAsMap();

        Assert.assertEquals(actual.size(), expectedReadContextData.size());
        for (KV<GATKRead, ReadContextData> kv : expectedReadContextData) {
            ReadContextData readContextData = actual.get(kv.getKey());
            Assert.assertNotNull(readContextData);
            Assert.assertTrue(CollectionUtils.isEqualCollection(Lists.newArrayList(readContextData.getOverlappingVariants()),
                    Lists.newArrayList(kv.getValue().getOverlappingVariants())));
            SimpleInterval minimalInterval = kv.getValue().getOverlappingReferenceBases().getInterval();
            ReferenceBases subset = readContextData.getOverlappingReferenceBases().getSubset(minimalInterval);
            Assert.assertEquals(subset, kv.getValue().getOverlappingReferenceBases());
        }
    }

    static class ReferenceBasesAnswer implements Answer<ReferenceBases>, Serializable {
        private static final long serialVersionUID = 1L;
        @Override
        public ReferenceBases answer(InvocationOnMock invocation) throws Throwable {
            return FakeReferenceSource.bases(invocation.getArgumentAt(1, SimpleInterval.class));
        }
    }
}