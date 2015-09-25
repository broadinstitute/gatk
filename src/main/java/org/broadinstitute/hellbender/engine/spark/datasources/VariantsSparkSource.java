package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.LongWritable;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;
import org.seqdoop.hadoop_bam.VCFInputFormat;
import org.seqdoop.hadoop_bam.VariantContextWritable;

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
     * Loads variants in parallel using Hadoop-BAM works for vcfs and bcfs.
     * @param vcfs files to load variants from.
     * @return JavaRDD<Variant> of variants from all files.
     */
    public JavaRDD<Variant> getParallelVariants(final List<String> vcfs) {
        throw new GATKException("This method does not currently work (issue with union()");
        /*
        JavaRDD<Variant> rddVariants = ctx.emptyRDD();
        for (String vcf : vcfs) {
            JavaRDD<Variant> variants = getParallelVariants(vcf);
            rddVariants.union(variants);
        }
        return rddVariants;
        */
    }

    /**
     * Loads variants in parallel using Hadoop-BAM for vcfs and bcfs.
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

}