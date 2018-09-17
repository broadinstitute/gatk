package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.collections.CollectionUtils;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
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
import java.util.function.Function;

public final class VariantsSparkSourceUnitTest extends GATKBaseTest {
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
                variantsSparkSource.getParallelVariants(vcf, null);

        List<GATKVariant> serialVariants = getSerialVariants(vcf);
        List<GATKVariant> parallelVariants = rddParallelVariants.collect();
        Assert.assertEquals(parallelVariants, serialVariants);
    }

    @Test(dataProvider = "loadVariants", groups = "spark")
    public void pairReadsAndVariantsTest_variantContexts(String vcf) {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        VariantsSparkSource variantsSparkSource = new VariantsSparkSource(ctx);
        JavaRDD<VariantContext> rddParallelVariantContexts =
                variantsSparkSource.getParallelVariantContexts(vcf, null);

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
                variantsSparkSource.getParallelVariants(vcfList, null);

        // retrieve the same set of variants, but through VariantsSource, and wrapped by
        // the same wrapper class used by VariantsSparkSource to facilitate comparison
        List<GATKVariant> variantsList =
                getVariantsListAs(vcfList, vc -> VariantContextVariantAdapter.sparkVariantAdapter(vc));

        Assert.assertTrue(CollectionUtils.isEqualCollection(rddParallelVariants.collect(), variantsList));
    }

    /**
     * getVariantsListAs grabs the variants from local files (or perhaps eventually buckets), applies
     * the wrapper function to each object, and returns them as a list of objects of the type returned
     * by the wrapper function
     * @param variantSources list of files  to read from
     * @param wrapFunction function applied to each VariantContext returned
     */
    private static <T> List<T> getVariantsListAs( List<String> variantSources, Function<VariantContext, T> wrapFunction ) {
        final List<T> aggregatedResults = new ArrayList<>();

        for ( final String variantSource : variantSources ) {
            try ( final FeatureDataSource<VariantContext> dataSource =
                          new FeatureDataSource<>(variantSource, null, 0, VariantContext.class) ) {
                aggregatedResults.addAll(wrapQueryResults(dataSource.iterator(), wrapFunction));
            }
        }
        return aggregatedResults;
    }

    private static <T> List<T> wrapQueryResults( final Iterator<VariantContext> queryResults, final Function<VariantContext, T> wrapFunction) {
        final List<T> wrappedResults = new ArrayList<>();
        while ( queryResults.hasNext() ) {
            wrappedResults.add(wrapFunction.apply(queryResults.next()));
        }
        return wrappedResults;
    }

    /**
     * Loads variants using FeatureDataSource<VariantContext>.
     * @param vcf file to load
     * @return List of GATKVariants.
     */
    static List<GATKVariant> getSerialVariants(final String vcf) {
        try ( final FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<>(new File(vcf), null, 0) ) {
            return Lists.newArrayList(wrapQueryResults(dataSource.iterator()));
        }
    }

    /**
     * Loads variants using FeatureDataSource<VariantContext>.
     * @param vcf file to load
     * @return List of VariantContext.
     */
    static List<VariantContext> getSerialVariantContexts(final String vcf) {
        try ( final FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<>(new File(vcf), null, 0) ) {
            return Lists.newArrayList(dataSource.iterator());
        }
    }

    static List<GATKVariant> wrapQueryResults(final Iterator<VariantContext> queryResults ) {
        final List<GATKVariant> wrappedResults = new ArrayList<>();
        while ( queryResults.hasNext() ) {
            wrappedResults.add(VariantContextVariantAdapter.sparkVariantAdapter(queryResults.next()));
        }
        return wrappedResults;
    }
}
