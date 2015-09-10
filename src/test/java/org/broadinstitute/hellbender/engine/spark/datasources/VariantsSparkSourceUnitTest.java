package org.broadinstitute.hellbender.engine.spark.datasources;

import com.beust.jcommander.internal.Lists;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class VariantsSparkSourceUnitTest extends BaseTest {
    @DataProvider(name = "loadVariants")
    public Object[][] loadVariants() {
        return new Object[][]{
                {hg19_chr1_1M_dbSNP},
                {hg19_chr1_1M_dbSNP_modified},
        };
    }

    @Test(dataProvider = "loadVariants", groups = "spark")
    public void pairReadsAndVariantsTest(String vcf) {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        VariantsSparkSource variantsSparkSource = new VariantsSparkSource(ctx);
        JavaRDD<Variant> rddSerialVariants =
                getSerialVariants(ctx, vcf);
        JavaRDD<Variant> rddParallelVariants =
                variantsSparkSource.getParallelVariants(vcf);

        List<Variant> serialVariants = rddSerialVariants.collect();
        List<Variant> parallelVariants = rddParallelVariants.collect();
        Assert.assertEquals(parallelVariants, serialVariants);
    }

    /**
     * Loads variants using FeatureDataSource<VariantContext>.
     * @param vcf file to load
     * @return JavaRDD of Variants.
     */
    static JavaRDD<Variant> getSerialVariants(final JavaSparkContext ctx, final String vcf) {
        List<Variant> records = Lists.newArrayList();
        try ( final FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<>(new File(vcf), getCodecForVariantSource(vcf), null, 0) ) {
            records.addAll(wrapQueryResults(dataSource.iterator()));
        }
        return ctx.parallelize(records);
    }

    @SuppressWarnings("unchecked")
    static FeatureCodec<VariantContext, ?> getCodecForVariantSource( final String variantSource ) {
        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(new File(variantSource));
        if ( !VariantContext.class.isAssignableFrom(codec.getFeatureType()) ) {
            throw new UserException(variantSource + " is not in a format that produces VariantContexts");
        }
        return (FeatureCodec<VariantContext, ?>)codec;
    }

    static List<Variant> wrapQueryResults( final Iterator<VariantContext> queryResults ) {
        final List<Variant> wrappedResults = new ArrayList<>();
        while ( queryResults.hasNext() ) {
            wrappedResults.add(VariantContextVariantAdapter.sparkVariantAdapter(queryResults.next()));
        }
        return wrappedResults;
    }
}