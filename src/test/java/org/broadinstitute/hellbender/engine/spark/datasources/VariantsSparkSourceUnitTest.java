package org.broadinstitute.hellbender.engine.spark.datasources;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkTestUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;
import java.util.UUID;

import static org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapterTest.createVariantContextVariantAdapterForTesting;

public class VariantsSparkSourceUnitTest extends BaseTest {
    @DataProvider(name = "loadVariants")
    public Object[][] loadVariants() {
        return new Object[][]{
                {hg19_chr1_1M_dbSNP},
                {hg19_chr1_1M_dbSNP_modified},
        };
    }

    @Test(dataProvider = "loadVariants")
    public void pairReadsAndVariantsTest(String vcf) {
        JavaSparkContext ctx = SparkTestUtils.getTestContext();

        VariantsSparkSource variantsSparkSource = new VariantsSparkSource(ctx);
        JavaRDD<Variant> rddSerialVariants =
                variantsSparkSource.getSerialVariants(vcf);
        JavaRDD<Variant> rddParallelVariants =
                variantsSparkSource.getParallelVariants(vcf);

        List<Variant> serialVariants = rddSerialVariants.collect();
        List<Variant> parallelVariants = rddParallelVariants.collect();
        Assert.assertEquals(parallelVariants, serialVariants);

        ctx.stop();
    }

}