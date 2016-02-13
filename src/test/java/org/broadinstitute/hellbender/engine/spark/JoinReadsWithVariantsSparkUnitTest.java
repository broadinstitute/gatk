package org.broadinstitute.hellbender.engine.spark;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMRecord;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class JoinReadsWithVariantsSparkUnitTest extends BaseTest {
    @DataProvider(name = "pairedReadsAndVariants")
    public Object[][] pairedReadsAndVariants(){
        List<Object[]> testCases = new ArrayList<>();

        for ( Class<?> readImplementation : Arrays.asList(Read.class, SAMRecord.class) ) {
            ReadsPreprocessingPipelineSparkTestData testData = new ReadsPreprocessingPipelineSparkTestData(readImplementation);
            List<GATKRead> reads = testData.getReads();
            List<Variant> variantList = testData.getVariants();
            List<KV<GATKRead, Iterable<Variant>>> kvReadiVariant = testData.getKvReadiVariant();

            testCases.add(new Object[]{reads, variantList, kvReadiVariant});
        }
        return testCases.toArray(new Object[][]{});
    }

    @Test(dataProvider = "pairedReadsAndVariants", groups = "spark")
    public void pairReadsAndVariantsTest(List<GATKRead> reads, List<Variant> variantList, List<KV<GATKRead, Iterable<Variant>>> kvReadiVariant) {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        JavaRDD<GATKRead> rddReads = ctx.parallelize(reads);
        JavaRDD<Variant> rddVariants = ctx.parallelize(variantList);
        JavaPairRDD<GATKRead, Iterable<Variant>> actual = BroadcastJoinReadsWithVariants.join(rddReads, rddVariants);
        Map<GATKRead, Iterable<Variant>> gatkReadIterableMap = actual.collectAsMap();

        Assert.assertEquals(gatkReadIterableMap.size(), kvReadiVariant.size());
        for (KV<GATKRead, Iterable<Variant>> kv : kvReadiVariant) {
            List<Variant> variants = Lists.newArrayList(gatkReadIterableMap.get(kv.getKey()));
            Assert.assertTrue(variants.stream().noneMatch( v -> v == null));
            HashSet<Variant> hashVariants = new HashSet<>(variants);
            final Iterable<Variant> iVariants = kv.getValue();
            HashSet<Variant> expectedHashVariants = Sets.newHashSet(iVariants);
            Assert.assertEquals(hashVariants, expectedHashVariants);
        }
    }
}