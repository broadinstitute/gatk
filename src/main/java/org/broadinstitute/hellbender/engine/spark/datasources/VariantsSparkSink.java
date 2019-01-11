package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFBlockCombiner;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFBlockCombiningIterator;
import org.disq_bio.disq.HtsjdkVariantsRdd;
import org.disq_bio.disq.HtsjdkVariantsRddStorage;
import scala.Tuple2;

import java.io.IOException;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

/**
 * VariantsSparkSink writes variants to a VCF file in parallel using Hadoop-BAM. BCF is not supported.
 */
public final class VariantsSparkSink {
    /**
     * Write variants to the given output file in VCF format with the given header. Note that writing sharded output is not supported.
     * @param ctx the JavaSparkContext
     * @param outputFile path to the output VCF
     * @param variants variants to write
     * @param header the header to put at the top of the output file
     * @throws IOException if an error occurs while writing
     */
    public static void writeVariants(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<VariantContext> variants,
            final VCFHeader header) throws IOException {
        writeVariants(ctx, outputFile, variants, header, false, null, 0, 0);
    }

    public static void writeVariants(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<VariantContext> variants,
            final VCFHeader header, final boolean writeGvcf, final List<Number> gqPartitions, final int defaultPloidy) throws IOException {
        writeVariants(ctx, outputFile, variants, header, writeGvcf, gqPartitions, defaultPloidy, 0);
    }

    /**
     * Write variants to the given output file in VCF format with the given header. Note that writing sharded output is not supported.
     * @param ctx the JavaSparkContext
     * @param outputFile path to the output VCF
     * @param variants variants to write
     * @param header the header to put at the top of the output file
     * @param numReducers the number of reducers to use when writing a single file. A value of zero indicates that the default
     *                    should be used.
     * @throws IOException if an error occurs while writing
     */
    public static void writeVariants(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<VariantContext> variants,
            final VCFHeader header, final boolean writeGvcf, final List<Number> gqPartitions, final int defaultPloidy,
            final int numReducers) throws IOException {
        String absoluteOutputFile = BucketUtils.makeFilePathAbsolute(outputFile);
        writeVariantsSingle(ctx, absoluteOutputFile, variants, header, writeGvcf, gqPartitions, defaultPloidy, numReducers);
    }

    private static void writeVariantsSingle(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<VariantContext> variants,
            final VCFHeader header, final boolean writeGvcf, final List<Number> gqPartitions, final int defaultPloidy, final int numReducers) throws IOException {

        //TODO remove me when https://github.com/broadinstitute/gatk/issues/4303 is fixed
        if (outputFile.endsWith(IOUtil.BCF_FILE_EXTENSION) || outputFile.endsWith(IOUtil.BCF_FILE_EXTENSION + ".gz")) {
            throw new UserException.UnimplementedFeature("It is currently not possible to write a BCF file on spark.  See https://github.com/broadinstitute/gatk/issues/4303 for more details .");
        }
        final JavaRDD<VariantContext> sortedVariants = sortVariants(variants, header, numReducers);
        final JavaRDD<VariantContext> variantsToSave;
        if (writeGvcf) {
            GVCFBlockCombiner gvcfBlockCombiner = new GVCFBlockCombiner(gqPartitions, defaultPloidy);
            gvcfBlockCombiner.addRangesToHeader(header);
            variantsToSave = sortedVariants.mapPartitions((FlatMapFunction<Iterator<VariantContext>, VariantContext>) v -> new GVCFBlockCombiningIterator(v, gqPartitions, defaultPloidy));
        } else {
            variantsToSave = sortedVariants;
        }
        HtsjdkVariantsRdd htsjdkVariantsRdd = new HtsjdkVariantsRdd(header, variantsToSave);
        HtsjdkVariantsRddStorage.makeDefault(ctx)
                .write(htsjdkVariantsRdd, outputFile);
    }

    private static JavaRDD<VariantContext> sortVariants(final JavaRDD<VariantContext> variants, final VCFHeader header, final int numReducers) {
        // Turn into key-value pairs so we can sort (by key). Values are null so there is no overhead in the amount
        // of data going through the shuffle.
        final JavaPairRDD<VariantContext, Void> rddVariantPairs = variants.mapToPair(variant -> new Tuple2<>(variant, (Void) null));

        // do a total sort so that all the records in partition i are less than those in partition i+1
        final Comparator<VariantContext> comparator = header.getVCFRecordComparator();
        final JavaPairRDD<VariantContext, Void> variantVoidPairs;
        if (comparator == null){
            variantVoidPairs = rddVariantPairs; //no sort
        } else if (numReducers > 0) {
            variantVoidPairs = rddVariantPairs.sortByKey(comparator, true, numReducers);
        } else {
            variantVoidPairs = rddVariantPairs.sortByKey(comparator);
        }

        return variantVoidPairs.map(Tuple2::_1);
    }
}
