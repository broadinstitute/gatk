package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.function.FlatMapFunction2;
import org.apache.spark.api.java.function.Function2;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFBlockCombiner;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFBlockCombiningIterator;
import org.disq_bio.disq.HtsjdkVariantsRdd;
import org.disq_bio.disq.HtsjdkVariantsRddStorage;
import scala.Tuple2;

import java.io.IOException;
import java.util.*;

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
            JavaRDD<VariantContext> sortedVariantsGvcf = sortedVariants.mapPartitions((FlatMapFunction<Iterator<VariantContext>, VariantContext>) v -> new GVCFBlockCombiningIterator(v, gqPartitions, defaultPloidy));
            variantsToSave = fixGvcfBlockBoundaries(ctx, sortedVariantsGvcf);
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

    /**
     * Remove GVCF block boundary artifacts that were introduced due to variants being processed in partitions.
     * The idea is that a hom-ref variant with an END attribute that is the last variant in a partition may be combined
     * with the first variant in the next partition in some circumstances (same GVCF band, both hom-ref). This method
     * makes the first variant in each partition available to the previous partition and then uses GVCFBlockCombiner
     * to combine the variant blocks if possible.
     */
    private static JavaRDD<VariantContext> fixGvcfBlockBoundaries(final JavaSparkContext ctx, final JavaRDD<VariantContext> variants) {
        return zipWithFirstElementInNextPartition(ctx, variants, new FlatMapFunction2<Iterator<VariantContext>, IndexedValue<VariantContext>, VariantContext>() {
            private static final long serialVersionUID = 1L;
            @Override
            public Iterator<VariantContext> call(Iterator<VariantContext> it, IndexedValue<VariantContext> indexedValue) {
                if (!it.hasNext()) {
                    return Collections.emptyIterator();
                }
                int startIndex = (indexedValue.index() == 1) ? 0 : 1; // first partition keeps first elt, others do not (note index is partition number elt comes from)
                // TODO: make this code more transparent (use more collections abstractions)
                // TODO: don't materialize iterator as list
                List<VariantContext> elts = Lists.newArrayList(it);
                elts = elts.subList(startIndex, elts.size());
                if (indexedValue.isPresent()) {
                    List<VariantContext> lastVcs = GVCFBlockCombiner.combine(elts.get(elts.size() - 1), indexedValue.get());
                    return Iterators.concat(elts.subList(0, elts.size() - 1).iterator(), lastVcs.iterator());
                }
                return elts.iterator();
            }
        }).map(vc -> {
            // clear DPs field
            vc.getGenotypes().forEach(g -> g.getExtendedAttributes().remove("DPs"));
            return vc;
        }
        );
    }

    static class IndexedValue<T> {

        private int index;
        private final T value;

        public IndexedValue(int index, T value) {
            this.index = index;
            this.value = value;
        }

        public boolean isPresent() {
            return value != null;
        }

        public T get() {
            return value;
        }

        public int index() {
            return index;
        }
    }

    private static <T> List<IndexedValue<T>> getFirstElementInEachPartition(final JavaRDD<T> rdd) {
        return rdd.mapPartitionsWithIndex((Function2<Integer, Iterator<T>, Iterator<IndexedValue<T>>>)
                (index, it) -> ImmutableList.of(new IndexedValue<>(index, it.hasNext() ? it.next() : null)).iterator(), true).collect();
    }

    /**
     * Allow access to the first element in the following partition of an RDD.
     * @param ctx the Spark context
     * @param rdd the RDD
     * @param fn a user function that processes the values in a partition while at the same time having access to the first element in the next partition
     * @param <T> the type of objects in the RDD
     * @return an RDD after the user function has been applied
     */
    private static <T> JavaRDD<T> zipWithFirstElementInNextPartition(final JavaSparkContext ctx, final JavaRDD<T> rdd, final FlatMapFunction2<Iterator<T>, IndexedValue<T>, T> fn) {
        int numPartitions = rdd.getNumPartitions();
        // index is partition it came *from*
        List<IndexedValue<T>> firstElementInEachPartition = getFirstElementInEachPartition(rdd);
        List<IndexedValue<T>> firstElementInEachPartitionShifted = new ArrayList<>(firstElementInEachPartition.subList(1, firstElementInEachPartition.size()));
        firstElementInEachPartitionShifted.add(new IndexedValue<>(numPartitions, null)); // next element for last partition is always null, since there are no more elements
        JavaRDD<IndexedValue<T>> firstElementInEachPartitionShiftedRdd = ctx.parallelize(firstElementInEachPartitionShifted, numPartitions);
        return rdd.zipPartitions(firstElementInEachPartitionShiftedRdd, (FlatMapFunction2<Iterator<T>, Iterator<IndexedValue<T>>, T>)
                (tIterator, indexedValueIterator) -> fn.call(tIterator, Iterators.getOnlyElement(indexedValueIterator)));
    }
}
