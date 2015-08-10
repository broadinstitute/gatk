package org.broadinstitute.hellbender.engine.spark.datasources;

import com.beust.jcommander.internal.Lists;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.LongWritable;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;
import org.seqdoop.hadoop_bam.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * VariantsSparkSource loads Variants from files serially (using FeatureDataSource<VariantContext>) or in parallel
 * using Hadoop-BAM.
 */
public final class VariantsSparkSource {
    private final JavaSparkContext ctx;

    public VariantsSparkSource(JavaSparkContext ctx) {
        this.ctx = ctx;
    }

    /**
     * Loads variants using FeatureDataSource<VariantContext>.
     * @param vcf file to load
     * @return JavaRDD of Variants.
     */
    public JavaRDD<Variant> getSerialVariants(final String vcf) {
        List<Variant> records = Lists.newArrayList();
        try ( final FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<>(new File(vcf), getCodecForVariantSource(vcf), null, 0) ) {
            records.addAll(wrapQueryResults(dataSource.iterator()));
        }
        return ctx.parallelize(records);
    }

    /**
     * Loads variants in parallel using Hadoop-BAM
     * @param vcfs files to load variants from.
     * @return JavaRDD<Variant> of variants from all files.
     */
    public JavaRDD<Variant> getParallelVariants(final List<String> vcfs) {
        JavaRDD<Variant> rddVariants = ctx.emptyRDD();
        for (String vcf : vcfs) {
            JavaRDD<Variant> variants = getParallelVariants(vcf);
            rddVariants.union(variants);
        }
        return rddVariants;
    }

    /**
     * Loads variants in parallel using Hadoop-BAM
     * @param vcf file to load variants from.
     * @return JavaRDD<Variant> of variants from all files.
     */
    public JavaRDD<Variant> getParallelVariants(final String vcf) {
        JavaPairRDD<LongWritable, VariantContextWritable> rdd2 = ctx.newAPIHadoopFile(
                vcf, VCFInputFormat.class, LongWritable.class, VariantContextWritable.class,
                new Configuration());
        return rdd2.map(v1 -> {
            VariantContext variant = v1._2().get();
            if (variant.getCommonInfo() == null) {
                throw new GATKException("no common info");
            }
            return VariantContextVariantAdapter.sparkVariantAdapter(variant);
        });
    }

    @SuppressWarnings("unchecked")
    private FeatureCodec<VariantContext, ?> getCodecForVariantSource( final String variantSource ) {
        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(new File(variantSource));
        if ( !VariantContext.class.isAssignableFrom(codec.getFeatureType()) ) {
            throw new UserException(variantSource + " is not in a format that produces VariantContexts");
        }
        return (FeatureCodec<VariantContext, ?>)codec;
    }

    private List<Variant> wrapQueryResults( final Iterator<VariantContext> queryResults ) {
        final List<Variant> wrappedResults = new ArrayList<>();
        while ( queryResults.hasNext() ) {
            wrappedResults.add(VariantContextVariantAdapter.sparkVariantAdapter(queryResults.next()));
        }
        return wrappedResults;
    }
}