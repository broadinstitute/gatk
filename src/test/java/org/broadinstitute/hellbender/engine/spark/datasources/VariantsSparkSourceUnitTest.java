package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.common.collect.Lists;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
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
import java.util.stream.Collectors;

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
        JavaRDD<Variant> rddParallelVariants =
                variantsSparkSource.getParallelVariants(vcf);

        List<Variant> serialVariants = getSerialVariants(vcf);
        List<Variant> parallelVariants = rddParallelVariants.collect();
        Assert.assertEquals(parallelVariants, serialVariants);
    }

    @Test(dataProvider = "loadVariants", groups = "spark")
    public void pairReadsAndVariantsTest_variantContexts(String vcf) {
        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        VariantsSparkSource variantsSparkSource = new VariantsSparkSource(ctx);
        JavaRDD<VariantContext> rddParallelVariants =
                variantsSparkSource.getParallelVariantContexts(vcf);

        List<VariantContext> serialVariants = getSerialVariantContexts(vcf);
        List<VariantContext> parallelVariants = rddParallelVariants.collect();

        final List<String> contigs = serialVariants.stream().map(vc -> vc.getContig()).distinct().collect(Collectors.toList());
        assertEquals(parallelVariants, serialVariants, new VariantContextComparator(contigs));
    }

    private void assertEquals(List<VariantContext> v1, List<VariantContext> v2, VariantContextComparator comparator) {
        if (v1.size() != v2.size()){
            throw new AssertionError("different sizes " + v1.size()+ " vs " + v2.size());
        }
        for (int i = 0; i < v1.size(); i++) {
            if (0 != comparator.compare(v1.get(i), v2.get(i))){
                throw new AssertionError("different element " + i + " " + v1.get(i) + " vs " + v2.get(i));
            }
        }
    }

    /**
     * Loads variants using FeatureDataSource<VariantContext>.
     * @param vcf file to load
     * @return List of Variants.
     */
    static List<Variant> getSerialVariants(final String vcf) {
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

    static List<Variant> wrapQueryResults( final Iterator<VariantContext> queryResults ) {
        final List<Variant> wrappedResults = new ArrayList<>();
        while ( queryResults.hasNext() ) {
            wrappedResults.add(VariantContextVariantAdapter.sparkVariantAdapter(queryResults.next()));
        }
        return wrappedResults;
    }
}