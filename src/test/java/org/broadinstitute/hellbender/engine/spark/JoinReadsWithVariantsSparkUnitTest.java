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
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class JoinReadsWithVariantsSparkUnitTest extends BaseTest {
    @DataProvider(name = "pairedReadsAndVariants")
    public Object[][] pairedReadsAndVariants(){
        List<Object[]> testCases = new ArrayList<>();

        for ( JoinStrategy joinStrategy : JoinStrategy.values() ) {
            for ( Class<?> readImplementation : Arrays.asList(Read.class, SAMRecord.class) ) {
                ReadsPreprocessingPipelineSparkTestData testData = new ReadsPreprocessingPipelineSparkTestData(readImplementation);
                List<GATKRead> reads = testData.getReads();
                List<GATKVariant> variantList = testData.getVariants();
                List<KV<GATKRead, Iterable<GATKVariant>>> kvReadiVariant = testData.getKvReadiVariant();

                testCases.add(new Object[]{reads, variantList, kvReadiVariant, joinStrategy});
            }
        }
        return testCases.toArray(new Object[][]{});
    }

    @Test(dataProvider = "pairedReadsAndVariants", groups = "spark")
    public void pairReadsAndVariantsTest(List<GATKRead> reads, List<GATKVariant> variantList, List<KV<GATKRead, Iterable<GATKVariant>>> kvReadiVariant, JoinStrategy joinStrategy) {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        JavaRDD<GATKRead> rddReads = ctx.parallelize(reads);
        JavaRDD<GATKVariant> rddVariants = ctx.parallelize(variantList);
        JavaPairRDD<GATKRead, Iterable<GATKVariant>> actual = joinStrategy == JoinStrategy.SHUFFLE ?
                                                          ShuffleJoinReadsWithVariants.join(rddReads, rddVariants) :
                                                          BroadcastJoinReadsWithVariants.join(rddReads, rddVariants, null);
        Map<GATKRead, Iterable<GATKVariant>> gatkReadIterableMap = actual.collectAsMap();

        Assert.assertEquals(gatkReadIterableMap.size(), kvReadiVariant.size());
        for (KV<GATKRead, Iterable<GATKVariant>> kv : kvReadiVariant) {
            List<GATKVariant> variants = Lists.newArrayList(gatkReadIterableMap.get(kv.getKey()));
            Assert.assertTrue(variants.stream().noneMatch( v -> v == null));
            HashSet<GATKVariant> hashVariants = new LinkedHashSet<>(variants);
            final Iterable<GATKVariant> iVariants = kv.getValue();
            HashSet<GATKVariant> expectedHashVariants = Sets.newLinkedHashSet(iVariants);
            Assert.assertEquals(hashVariants, expectedHashVariants);
        }
    }
}