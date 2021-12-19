package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.util.FileExtensions;
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
import org.disq_bio.disq.TabixIndexWriteOption;
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
        writeVariants(ctx, outputFile, variants, header, false, null, 0, 0, true, true);
    }

    /**
     * Write variants to the given output file in VCF format with the given header. Note that writing sharded output is not supported.
     * @param ctx the JavaSparkContext
     * @param outputFile path to the output VCF
     * @param variants variants to write
     * @param header the header to put at the top of the output file
     * @param writeTabixIndex whether to write a tabix index (bgzipped VCF only)
     * @throws IOException if an error occurs while writing
     */
    public static void writeVariants(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<VariantContext> variants,
            final VCFHeader header, final boolean writeTabixIndex) throws IOException {
        writeVariants(ctx, outputFile, variants, header, writeTabixIndex, true);
    }

    /**
     * Write variants to the given output file in VCF format with the given header. Note that writing sharded output is not supported.
     * @param ctx the JavaSparkContext
     * @param outputFile path to the output VCF
     * @param variants variants to write
     * @param header the header to put at the top of the output file
     * @param writeTabixIndex whether to write a tabix index (bgzipped VCF only)
     * @param sortVariantsToHeader whether to sort variants to be consistent with the header
     * @throws IOException if an error occurs while writing
     */
    public static void writeVariants(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<VariantContext> variants,
            final VCFHeader header, final boolean writeTabixIndex, final boolean sortVariantsToHeader) throws IOException {
        writeVariants(ctx, outputFile, variants, header, false, null, 0, 0, writeTabixIndex, sortVariantsToHeader);
    }

    /**
     * Write variants to the given output file in VCF format with the given header. Note that writing sharded output is not supported.
     * @param ctx the JavaSparkContext
     * @param outputFile path to the output VCF
     * @param variants variants to write
     * @param header the header to put at the top of the output file
     * @param writeGvcf whether to write GVCF output
     * @param gqPartitions the GQ partitions for GVCF output
     * @param defaultPloidy the default ploidy for GVCF output
     * @throws IOException if an error occurs while writing
     */
    public static void writeVariants(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<VariantContext> variants,
            final VCFHeader header, final boolean writeGvcf, final List<Number> gqPartitions, final int defaultPloidy) throws IOException {
        writeVariants(ctx, outputFile, variants, header, writeGvcf, gqPartitions, defaultPloidy, 0, true, true);
    }

    /**
     * Write variants to the given output file in VCF format with the given header. Note that writing sharded output is not supported.
     * @param ctx the JavaSparkContext
     * @param outputFile path to the output VCF
     * @param variants variants to write
     * @param header the header to put at the top of the output file
     * @param writeGvcf whether to write GVCF output
     * @param gqPartitions the GQ partitions for GVCF output
     * @param defaultPloidy the default ploidy for GVCF output
     * @param numReducers the number of reducers to use when writing a single file. A value of zero indicates that the default
     *                    should be used.
     * @param writeTabixIndex whether to write a tabix index (bgzipped VCF only)
     * @param sortVariantsToHeader whether to sort variants to be consistent with the header
     * @throws IOException if an error occurs while writing
     */
    public static void writeVariants(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<VariantContext> variants,
            final VCFHeader header, final boolean writeGvcf, final List<Number> gqPartitions, final int defaultPloidy,
            final int numReducers, final boolean writeTabixIndex, final boolean sortVariantsToHeader) throws IOException {
        String absoluteOutputFile = BucketUtils.makeFilePathAbsolute(outputFile);
        writeVariantsSingle(ctx, absoluteOutputFile, variants, header, writeGvcf, gqPartitions, defaultPloidy, numReducers, writeTabixIndex, sortVariantsToHeader);
    }

    private static void writeVariantsSingle(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<VariantContext> variants,
            final VCFHeader header, final boolean writeGvcf, final List<Number> gqPartitions, final int defaultPloidy,
            final int numReducers, final boolean writeTabixIndex, final boolean sortVariantsToHeader) throws IOException {

        //TODO remove me when https://github.com/broadinstitute/gatk/issues/4303 is fixed
        if (outputFile.endsWith(FileExtensions.BCF) || outputFile.endsWith(FileExtensions.BCF + ".gz")) {
            throw new UserException.UnimplementedFeature("It is currently not possible to write a BCF file on spark.  See https://github.com/broadinstitute/gatk/issues/4303 for more details .");
        }
        final JavaRDD<VariantContext> sortedVariants = sortVariantsToHeader ? sortVariants(variants, header, numReducers) : variants;
        final JavaRDD<VariantContext> variantsToSave;
        if (writeGvcf) {
            GVCFBlockCombiner gvcfBlockCombiner = new GVCFBlockCombiner(gqPartitions, false);
            gvcfBlockCombiner.addRangesToHeader(header);
            variantsToSave = sortedVariants.mapPartitions((FlatMapFunction<Iterator<VariantContext>, VariantContext>) v -> new GVCFBlockCombiningIterator(v, gqPartitions, defaultPloidy));
        } else {
            variantsToSave = sortedVariants;
        }
        TabixIndexWriteOption tabixIndexWriteOption = TabixIndexWriteOption.fromBoolean(writeTabixIndex);
        HtsjdkVariantsRdd htsjdkVariantsRdd = new HtsjdkVariantsRdd(header, variantsToSave);
        HtsjdkVariantsRddStorage.makeDefault(ctx)
                .write(htsjdkVariantsRdd, outputFile, tabixIndexWriteOption);
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
