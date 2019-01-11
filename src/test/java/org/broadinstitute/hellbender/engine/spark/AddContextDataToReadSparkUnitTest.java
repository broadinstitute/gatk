package org.broadinstitute.hellbender.engine.spark;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReadContextData;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.utils.KV;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.testutils.FakeReferenceSource;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.Serializable;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class AddContextDataToReadSparkUnitTest extends GATKBaseTest {
    @DataProvider(name = "bases")
    public Object[][] bases() {
        List<Class<?>> classes = Collections.singletonList(SAMRecord.class);
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
                                   JoinStrategy joinStrategy) {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        JavaRDD<GATKRead> rddReads = ctx.parallelize(reads);
        JavaRDD<GATKVariant> rddVariants = ctx.parallelize(variantList);

        SAMSequenceDictionary sd = new SAMSequenceDictionary(Lists.newArrayList(new SAMSequenceRecord("1", 100000), new SAMSequenceRecord("2", 100000)));
        JavaPairRDD<GATKRead, ReadContextData> rddActual = AddContextDataToReadSpark.add(ctx, rddReads,
                new TestMultiReferenceSource(sd),
                rddVariants, null, joinStrategy,
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

    // Provide a fake implementation of this class for testing. We can't use a real mock since this is used as a Spark
    // broadcast variable. Mocks are mutated when they're accessed, which can result in ConcurrentModificationExceptions
    // during serialization/broadcast.
    static class TestMultiReferenceSource extends ReferenceMultiSparkSource implements Serializable {
        private static final long serialVersionUID = 1L;

        final SAMSequenceDictionary sequenceDictionary;

        public TestMultiReferenceSource(final SAMSequenceDictionary sd) {
            sequenceDictionary = sd;
        }

        @Override
        public ReferenceBases getReferenceBases(final SimpleInterval interval) {
            return FakeReferenceSource.bases(interval);
        }

        @Override
        public SAMSequenceDictionary getReferenceSequenceDictionary(final SAMSequenceDictionary optReadSequenceDictionaryToMatch) {
            return sequenceDictionary;
        }

        @Override
        public boolean isCompatibleWithSparkBroadcast(){
            return true;
        }

        @Override
        public SerializableFunction<GATKRead, SimpleInterval> getReferenceWindowFunction() {
            return ReferenceWindowFunctions.IDENTITY_FUNCTION;
        }

    }

}