package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.LongWritable;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;
import org.seqdoop.hadoop_bam.VCFInputFormat;
import org.seqdoop.hadoop_bam.VariantContextWritable;

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
     * Loads variants in parallel using Hadoop-BAM for vcfs and bcfs.
     * @param vcf file to load variants from.
     * @return JavaRDD<GATKVariant> of variants from all files.
     */
    public JavaRDD<GATKVariant> getParallelVariants(final String vcf) {
        return getParallelVariantContexts(vcf).filter(vc -> vc.getCommonInfo() != null).map(vc -> VariantContextVariantAdapter.sparkVariantAdapter(vc));
    }

    /**
     * Loads variants in parallel using Hadoop-BAM for vcfs and bcfs.
     * @param vcf file to load variants from.
     * @return JavaRDD<VariantContext> of variants from all files.
     */
    public JavaRDD<VariantContext> getParallelVariantContexts(final String vcf) {
        final JavaPairRDD<LongWritable, VariantContextWritable> rdd2 = ctx.newAPIHadoopFile(
                vcf, VCFInputFormat.class, LongWritable.class, VariantContextWritable.class,
                new Configuration());
        return rdd2.map(v1 -> v1._2().get());
    }

}