package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.common.collect.Lists;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.collections.CollectionUtils;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.engine.datasources.VariantsSource;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public final class VariantsSparkSourceUnitTest extends BaseTest {
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
        JavaRDD<GATKVariant> rddParallelVariants =
                variantsSparkSource.getParallelVariants(vcf);

        List<GATKVariant> serialVariants = getSerialVariants(vcf);
        List<GATKVariant> parallelVariants = rddParallelVariants.collect();
        Assert.assertEquals(parallelVariants, serialVariants);
    }

    @Test(dataProvider = "loadVariants", groups = "spark")
    public void pairReadsAndVariantsTest_variantContexts(String vcf) {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        VariantsSparkSource variantsSparkSource = new VariantsSparkSource(ctx);
        JavaRDD<VariantContext> rddParallelVariantContexts =
                variantsSparkSource.getParallelVariantContexts(vcf);

        VariantContextTestUtils.assertEqualVariants(getSerialVariantContexts(vcf), rddParallelVariantContexts.collect());
    }

    @DataProvider(name = "loadMultipleVCFs")
    public Object[][] loadMultipleVCFs() {
        return new Object[][]{
                {Arrays.asList(hg19_chr1_1M_dbSNP, hg19_chr1_1M_dbSNP_modified)},
        };
    }

    @Test(dataProvider = "loadMultipleVCFs", groups = "spark")
    public void getMultipleParallelVCFsTest(List<String> vcfList) {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        VariantsSparkSource variantsSparkSource = new VariantsSparkSource(ctx);

        JavaRDD<GATKVariant> rddParallelVariants =
                variantsSparkSource.getParallelVariants(vcfList);

        // retrieve the same set of variants, but through VariantsSource, and wrapped by
        // the same wrapper class used by VariantsSparkSource to facilitate comparison
        List<GATKVariant> variantsList =
                VariantsSource.getVariantsListAs
                        (vcfList, vc -> VariantContextVariantAdapter.sparkVariantAdapter(vc));

        Assert.assertTrue(CollectionUtils.isEqualCollection(rddParallelVariants.collect(), variantsList));
    }

    /**
     * Loads variants using FeatureDataSource<VariantContext>.
     * @param vcf file to load
     * @return List of GATKVariants.
     */
    static List<GATKVariant> getSerialVariants(final String vcf) {
        try ( final FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<>(new File(vcf), getCodecForVariantSource(vcf), null, 0) ) {
            return Lists.newArrayList(wrapQueryResults(dataSource.iterator()));
        }
    }

    /**
     * Loads variants using FeatureDataSource<VariantContext>.
     * @param vcf file to load
     * @return List of VariantContext.
     */
    static List<VariantContext> getSerialVariantContexts(final String vcf) {
        try ( final FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<>(new File(vcf), getCodecForVariantSource(vcf), null, 0) ) {
            return Lists.newArrayList(dataSource.iterator());
        }
    }

    @SuppressWarnings("unchecked")
    static FeatureCodec<VariantContext, ?> getCodecForVariantSource( final String variantSource ) {
        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(new File(variantSource));
        if ( !VariantContext.class.isAssignableFrom(codec.getFeatureType()) ) {
            throw new UserException(variantSource + " is not in a format that produces VariantContexts");
        }
        return (FeatureCodec<VariantContext, ?>)codec;
    }

    static List<GATKVariant> wrapQueryResults(final Iterator<VariantContext> queryResults ) {
        final List<GATKVariant> wrappedResults = new ArrayList<>();
        while ( queryResults.hasNext() ) {
            wrappedResults.add(VariantContextVariantAdapter.sparkVariantAdapter(queryResults.next()));
        }
        return wrappedResults;
    }
}
