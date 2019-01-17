package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.compress.CompressionCodec;
import org.apache.hadoop.mapred.FileAlreadyExistsException;
import org.apache.hadoop.mapreduce.JobContext;
import org.apache.hadoop.mapreduce.OutputFormat;
import org.apache.hadoop.mapreduce.RecordWriter;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFWriter;
import org.seqdoop.hadoop_bam.*;
import org.seqdoop.hadoop_bam.util.BGZFCodec;
import org.seqdoop.hadoop_bam.util.VCFFileMerger;
import scala.Tuple2;

import java.io.IOException;
import java.io.OutputStream;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * VariantsSparkSink writes variants to a VCF file in parallel using Hadoop-BAM. BCF is not supported.
 */
public final class VariantsSparkSink {

    // Output format class for writing VCF files through saveAsNewAPIHadoopFile. Must be public.
    public static class SparkVCFOutputFormat extends KeyIgnoringVCFOutputFormat<NullWritable> {
        public static VCFHeader vcfHeader;

        public static void setVCFHeader(final VCFHeader header) {
            vcfHeader = header;
        }

        public SparkVCFOutputFormat() {
            super(VCFFormat.VCF);
        }

        @Override
        public RecordWriter<NullWritable, VariantContextWritable> getRecordWriter(TaskAttemptContext ctx) throws IOException {
            setHeader(vcfHeader);
            // don't add an extension, since FileOutputFormat will add a compression extension automatically (e.g. .bgz)
            return super.getRecordWriter(ctx);
        }

        @Override
        public void checkOutputSpecs(JobContext job) throws IOException {
            try {
                super.checkOutputSpecs(job);
            } catch (FileAlreadyExistsException e) {
                // delete existing files before overwriting them
                final Path outDir = getOutputPath(job);
                outDir.getFileSystem(job.getConfiguration()).delete(outDir, true);
            }
        }
    }

    public static class SparkHeaderlessVCFOutputFormat extends SparkVCFOutputFormat {
        static final String GVCF = "hadoopbam.vcf.gvcf";
        static final String GQ_PARTITIONS = "hadoopbam.vcf.gvcf.gq-partitions";
        static final String DEFAULT_PLOIDY = "hadoopbam.vcf.gvcf.default-ploidy";

        public static void writeGvcf(Configuration conf, List<Integer> gqPartitions, int defaultPloidy) {
            conf.setBoolean(GVCF, true);
            conf.set(GQ_PARTITIONS, Joiner.on(",").join(gqPartitions));
            conf.setInt(DEFAULT_PLOIDY, defaultPloidy);
        }

        public static void unsetGvcf(Configuration conf) {
            conf.unset(GVCF);
            conf.unset(GQ_PARTITIONS);
            conf.unset(DEFAULT_PLOIDY);
        }

        @Override
        public RecordWriter<NullWritable, VariantContextWritable> getRecordWriter(TaskAttemptContext ctx, Path out) throws IOException {
            if (this.header == null) {
                throw new IOException("Can't create a RecordWriter without the VCF header");
            }

            final boolean wh = ctx.getConfiguration().getBoolean(
                    WRITE_HEADER_PROPERTY, true);

            switch (format) {
                case BCF: return new KeyIgnoringBCFRecordWriter<>(out, header, wh, ctx);
                case VCF: return new GvcfKeyIgnoringVCFRecordWriter<>(out, header, wh, ctx);
                default: throw new IllegalStateException("Unrecognized variant format: " + format);
            }
        }

        @Override
        public RecordWriter<NullWritable, VariantContextWritable> getRecordWriter(TaskAttemptContext ctx, OutputStream outputStream) throws IOException {
            if (this.header == null) {
                throw new IOException("Can't create a RecordWriter without the VCF header");
            }

            final boolean wh = ctx.getConfiguration().getBoolean(
                    WRITE_HEADER_PROPERTY, true);

            switch (format) {
                case BCF: return new KeyIgnoringBCFRecordWriter<>(outputStream, header, wh, ctx);
                case VCF: return new GvcfKeyIgnoringVCFRecordWriter<>(outputStream, header, wh, ctx);
                default: throw new IllegalStateException("Unrecognized variant format: " + format);
            }
        }

        @Override
        public RecordWriter<NullWritable, VariantContextWritable> getRecordWriter(TaskAttemptContext ctx) throws IOException {
            ctx.getConfiguration().setBoolean(WRITE_HEADER_PROPERTY, false);
            return super.getRecordWriter(ctx);
        }
    }

    public static class GvcfKeyIgnoringVCFRecordWriter<K> extends KeyIgnoringVCFRecordWriter<K> {

        public GvcfKeyIgnoringVCFRecordWriter(Path output, VCFHeader header, boolean writeHeader, TaskAttemptContext ctx) throws IOException {
            super(output, header, writeHeader, ctx);
        }

        public GvcfKeyIgnoringVCFRecordWriter(OutputStream outputStream, VCFHeader header, boolean writeHeader, TaskAttemptContext ctx) throws IOException {
            super(outputStream, header, writeHeader, ctx);
        }

        @Override
        protected VariantContextWriter createVariantContextWriter(Configuration conf,
                                                                  OutputStream out) {
            return getVariantContextWriter(super.createVariantContextWriter(conf, out), conf);
        }
    }

    /**
     * Wrap in a GVCF writer if GVCF output is configured.
     */
    private static VariantContextWriter getVariantContextWriter(VariantContextWriter writer, Configuration conf) {
        if (conf.getBoolean(SparkHeaderlessVCFOutputFormat.GVCF, false)) {
            List<Integer> gqPartitions = StreamSupport.stream(
                    Splitter.on(",").split(conf.get(SparkHeaderlessVCFOutputFormat.GQ_PARTITIONS)).spliterator(), false)
                    .map(Integer::parseInt)
                    .collect(Collectors.toList());
            int defaultPloidy = conf.getInt(SparkHeaderlessVCFOutputFormat.DEFAULT_PLOIDY, HomoSapiensConstants.DEFAULT_PLOIDY);
            return new GVCFWriter(writer, gqPartitions, defaultPloidy);
        }
        return writer;
    }

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
            final VCFHeader header, final boolean writeGvcf, final List<Integer> gqPartitions, final int defaultPloidy) throws IOException {
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
            final VCFHeader header, final boolean writeGvcf, final List<Integer> gqPartitions, final int defaultPloidy,
            final int numReducers) throws IOException {
        String absoluteOutputFile = BucketUtils.makeFilePathAbsolute(outputFile);
        writeVariantsSingle(ctx, absoluteOutputFile, variants, header, writeGvcf, gqPartitions, defaultPloidy, numReducers);
    }

    private static void writeVariantsSingle(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<VariantContext> variants,
            final VCFHeader header, final boolean writeGvcf, final List<Integer> gqPartitions, final int defaultPloidy, final int numReducers) throws IOException {

        final Configuration conf = ctx.hadoopConfiguration();

        //TODO remove me when https://github.com/broadinstitute/gatk/issues/4303 are fixed
        if (outputFile.endsWith(IOUtil.BCF_FILE_EXTENSION) || outputFile.endsWith(IOUtil.BCF_FILE_EXTENSION + ".gz")) {
            throw new UserException.UnimplementedFeature("It is currently not possible to write a BCF file on spark.  See https://github.com/broadinstitute/gatk/issues/4303 for more details .");
        }

        if (outputFile.endsWith(BGZFCodec.DEFAULT_EXTENSION) || outputFile.endsWith(".gz")) {
            conf.setBoolean(FileOutputFormat.COMPRESS, true);
            conf.setClass(FileOutputFormat.COMPRESS_CODEC, BGZFCodec.class, CompressionCodec.class);
        } else {
            conf.setBoolean(FileOutputFormat.COMPRESS, false);
        }
        if (writeGvcf) {
            SparkHeaderlessVCFOutputFormat.writeGvcf(conf, gqPartitions, defaultPloidy);
        } else {
            SparkHeaderlessVCFOutputFormat.unsetGvcf(conf);
        }

        final JavaRDD<VariantContext> sortedVariants = sortVariants(variants, header, numReducers);
        final String outputPartsDirectory = outputFile + ".parts/";
        saveAsShardedHadoopFiles(ctx, conf, outputPartsDirectory, sortedVariants,  header, false);
        VCFFileMerger.mergeParts(outputPartsDirectory, outputFile, header);
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

    private static void saveAsShardedHadoopFiles(
            final JavaSparkContext ctx, final Configuration conf, final String outputFile, JavaRDD<VariantContext> variants,
            final VCFHeader header, final boolean writeHeader) throws IOException {
        // Set the static header on the driver thread.
        SparkVCFOutputFormat.setVCFHeader(header);

        final Broadcast<VCFHeader> headerBroadcast = ctx.broadcast(header);

        // SparkVCFOutputFormat is a static class, so we need to copy the header to each worker then call
        final JavaRDD<VariantContext> variantsRDD = setHeaderForEachPartition(variants, headerBroadcast);

        // The expected format for writing is JavaPairRDD where the key is ignored and the value is VariantContextWritable.
        final JavaPairRDD<VariantContext, VariantContextWritable> rddVariantContextWriteable = pairVariantsWithVariantContextWritables(variantsRDD);

        rddVariantContextWriteable.saveAsNewAPIHadoopFile(outputFile, VariantContext.class, VariantContextWritable.class, getOutputFormat(writeHeader), conf);
    }

    private static JavaRDD<VariantContext> setHeaderForEachPartition(final JavaRDD<VariantContext> variants, final Broadcast<VCFHeader> headerBroadcast) {
        return variants.mapPartitions(iterator -> {
            SparkVCFOutputFormat.setVCFHeader(headerBroadcast.getValue());
            return iterator;
        });
    }

    private static JavaPairRDD<VariantContext, VariantContextWritable> pairVariantsWithVariantContextWritables(JavaRDD<VariantContext> records) {
        return records.mapToPair(variantContext -> {
            final VariantContextWritable variantContextWritable = new VariantContextWritable();
            variantContextWritable.set(variantContext);
            return new Tuple2<>(variantContext, variantContextWritable);
        });
    }

    private static Class<? extends OutputFormat<NullWritable, VariantContextWritable>> getOutputFormat(final boolean writeHeader) {
        return writeHeader ? SparkVCFOutputFormat.class : SparkHeaderlessVCFOutputFormat.class;
    }
}
