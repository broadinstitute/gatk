package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import org.apache.commons.collections4.iterators.IteratorIterable;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapred.FileAlreadyExistsException;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.JobContext;
import org.apache.hadoop.mapreduce.RecordWriter;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.parquet.avro.AvroParquetOutputFormat;
import org.apache.spark.RangePartitioner;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.bdgenomics.adam.models.RecordGroupDictionary;
import org.bdgenomics.adam.models.SequenceDictionary;
import org.bdgenomics.formats.avro.AlignmentRecord;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GATKReadToBDGAlignmentRecordConverter;
import org.broadinstitute.hellbender.utils.read.HeaderlessSAMRecordCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.seqdoop.hadoop_bam.KeyIgnoringBAMOutputFormat;
import org.seqdoop.hadoop_bam.SAMFormat;
import org.seqdoop.hadoop_bam.SAMRecordWritable;
import org.seqdoop.hadoop_bam.util.SAMOutputPreparer;
import scala.Tuple2;
import scala.math.Ordering;
import scala.math.Ordering$;
import scala.reflect.ClassTag;
import scala.reflect.ClassTag$;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Comparator;

/**
 * ReadsSparkSink writes GATKReads to a file. This code lifts from the HadoopGenomics/Hadoop-BAM
 * read writing code as well as from bigdatagenomics/adam.
 */
public final class ReadsSparkSink {
    private static final Logger logger = LogManager.getLogger(ReadsSparkSink.class);

    // We need an output format for saveAsNewAPIHadoopFile.
    public static class SparkBAMOutputFormat extends KeyIgnoringBAMOutputFormat<NullWritable> {
        public static SAMFileHeader bamHeader = null;

        public static void setHeader(final SAMFileHeader header) {
            bamHeader = header;
        }

        @Override
        public RecordWriter<NullWritable, SAMRecordWritable> getRecordWriter(TaskAttemptContext ctx) throws IOException {
            setSAMHeader(bamHeader);
            return super.getRecordWriter(ctx);
        }

        @Override
        public void checkOutputSpecs(JobContext job) throws IOException {
            try {
                super.checkOutputSpecs(job);
            } catch (FileAlreadyExistsException e) {
                // delete existing files before overwriting them
                Path outDir = getOutputPath(job);
                outDir.getFileSystem(job.getConfiguration()).delete(outDir, true);
            }
        }
    }

    public static class SparkHeaderlessBAMOutputFormat extends SparkBAMOutputFormat {
        public SparkHeaderlessBAMOutputFormat() {
            setWriteHeader(false);
        }
    }

    /**
     * writeReads writes rddReads to outputFile with header as the file header.
     * @param ctx the JavaSparkContext to write.
     * @param outputFile path to the output bam.
     * @param reads reads to write.
     * @param header the header to put at the top of the files
     * @param format should the output be a single file, sharded, ADAM, etc.
     */
    public static void writeReads(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<GATKRead> reads,
            final SAMFileHeader header, ReadsWriteFormat format) throws IOException {

        if (IOUtils.isCramFileName(outputFile)) {
            throw new UserException("Writing CRAM files in Spark tools is not supported yet.");
        }

        // The underlying reads are required to be in SAMRecord format in order to be
        // written out, so we convert them to SAMRecord explicitly here. If they're already
        // SAMRecords, this will effectively be a no-op. The SAMRecords will be headerless
        // for efficient serialization.
        final JavaRDD<SAMRecord> samReads = reads.map(read -> read.convertToSAMRecord(null));

        if (format.equals(ReadsWriteFormat.SINGLE)) {
            header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
            writeReadsSingle(ctx, outputFile, samReads, header);
        } else if (format.equals(ReadsWriteFormat.SHARDED)) {
            writeReadsSharded(ctx, outputFile, samReads, header);
        } else if (format.equals(ReadsWriteFormat.ADAM)) {
            writeReadsADAM(ctx, outputFile, samReads, header);
        }
    }

    private static void writeReadsADAM(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<SAMRecord> reads,
            final SAMFileHeader header) throws IOException {
        SequenceDictionary seqDict = SequenceDictionary.fromSAMSequenceDictionary(header.getSequenceDictionary());
        RecordGroupDictionary readGroups = RecordGroupDictionary.fromSAMHeader(header);
        JavaPairRDD<Void, AlignmentRecord> rddAlignmentRecords =
                reads.map(read -> {
                    read.setHeader(header);
                    AlignmentRecord alignmentRecord = GATKReadToBDGAlignmentRecordConverter.convert(read, seqDict, readGroups);
                    read.setHeader(null); // Restore the header to its previous state so as not to surprise the caller
                    return alignmentRecord;
                }).mapToPair(alignmentRecord -> new Tuple2<>(null, alignmentRecord));
        // instantiating a Job is necessary here in order to set the Hadoop Configuration...
        Job job = Job.getInstance(ctx.hadoopConfiguration());
        // ...here, which sets a config property that the AvroParquetOutputFormat needs when writing data. Specifically,
        // we are writing the Avro schema to the Configuration as a JSON string. The AvroParquetOutputFormat class knows
        // how to translate objects in the Avro data model to the Parquet primitives that get written.
        AvroParquetOutputFormat.setSchema(job, AlignmentRecord.getClassSchema());
        deleteHadoopFile(outputFile);
        rddAlignmentRecords.saveAsNewAPIHadoopFile(
                outputFile, Void.class, AlignmentRecord.class, AvroParquetOutputFormat.class, job.getConfiguration());
    }

    private static void writeReadsSharded(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<SAMRecord> reads,
            final SAMFileHeader header) throws IOException {
        // Set the header on the main thread.
        SparkBAMOutputFormat.setHeader(header);
        // MyOutputFormat is a static class, so we need to copy the header to each worker then call
        final Broadcast<SAMFileHeader> broadcast = ctx.broadcast(header);
        JavaRDD<SAMRecord> readsRDD = reads.mapPartitions(readIterator -> {
            SparkBAMOutputFormat.setHeader(broadcast.getValue());
            return new IteratorIterable<>(readIterator);
        });

        // The expected format for writing is JavaPairRDD where the key is ignored and the value is a
        // SAMRecordWritable. We use SAMRecord as the key, so we can sort the reads before writing them to a file.
        JavaPairRDD<SAMRecord, SAMRecordWritable> rddSamRecordWriteable = readsRDD.mapToPair(read -> {
            read.setHeader(broadcast.getValue());
            SAMRecordWritable samRecordWritable = new SAMRecordWritable();
            samRecordWritable.set(read);
            return new Tuple2<>(read, samRecordWritable);
        });

        deleteHadoopFile(outputFile);
        rddSamRecordWriteable.saveAsNewAPIHadoopFile(outputFile, SAMRecord.class, SAMRecordWritable.class, SparkBAMOutputFormat.class);
    }

    private static void writeReadsSingle(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<SAMRecord> reads,
            final SAMFileHeader header) throws IOException {
        // Set the header on the main thread.
        SparkBAMOutputFormat.setHeader(header);

        // Turn into key-value pairs so we can sort (by key). Values are null so there is no overhead in the amount
        // of data going through the shuffle.
        JavaPairRDD<SAMRecord, Void> rddReadPairs = reads.mapToPair(read -> new Tuple2<>(read, (Void) null));

        // do a total sort so that all the reads in partition i are less than those in partition i+1
        final JavaPairRDD<SAMRecord, Void> out =
                totalSortByKey(rddReadPairs, new HeaderlessSAMRecordCoordinateComparator(header), SAMRecord.class);

        // MyOutputFormat is a static class, so we need to copy the header to each worker then call
        // MyOutputFormat.setHeader.
        final Broadcast<SAMFileHeader> broadcast = ctx.broadcast(header);

        final JavaPairRDD<SAMRecord, SAMRecordWritable> finalOut = out.mapPartitions(tuple2Iterator -> {
            SparkBAMOutputFormat.setHeader(broadcast.getValue());
            return new IteratorIterable<>(tuple2Iterator);
        }).mapToPair(t -> {
            SAMRecord samRecord = t._1();
            samRecord.setHeader(header);
            SAMRecordWritable samRecordWritable = new SAMRecordWritable();
            samRecordWritable.set(samRecord);
            return new Tuple2<>(samRecord, samRecordWritable);
        });

        deleteHadoopFile(outputFile);
        // Use SparkHeaderlessBAMOutputFormat so that the header is not written for each part file
        finalOut.saveAsNewAPIHadoopFile(outputFile, SAMRecord.class, SAMRecordWritable.class, SparkHeaderlessBAMOutputFormat.class);
        mergeBam(outputFile, header);
    }

    /**
     * Perform a totally-ordered sort by the given RDD's keys, so that the keys within a partition are sorted and, in
     * addition, all the keys in partition <i>i</i> are less than those in partition <i>i+1</i>. The implementation uses
     * Spark's RangePartitioner to create roughly equal partitions by sampling the keys in the RDD.
     */
    public static <K, V> JavaPairRDD<K, V> totalSortByKey(JavaPairRDD<K, V> rdd, Comparator<K> comparator, Class<K> keyClass) {
        int numPartitions = rdd.partitions().size();
        if (numPartitions == 1) {
            return rdd.sortByKey(comparator);
        }
        logger.info("****NUM PARTITIONS: " + numPartitions);
        Ordering<K> ordering = Ordering$.MODULE$.comparatorToOrdering(comparator);
        ClassTag<K> tag = ClassTag$.MODULE$.apply(keyClass);
        RangePartitioner<K, V> partitioner = new RangePartitioner<>(numPartitions, rdd.rdd(), true, ordering, tag);
        return rdd.repartitionAndSortWithinPartitions(partitioner, comparator);
    }

    private static void deleteHadoopFile(String fileToObliterate) throws IOException {
        Configuration conf = new Configuration();
        FileSystem fs = FileSystem.get(conf);
        fs.delete(new Path(fileToObliterate), true);
    }

    private static void mergeBam(String outputFile, SAMFileHeader header) throws IOException {
        // At this point, the part files (part-r-00000, part-r-00001, etc) are in a directory named outputFile.
        // Each part file is a BAM file with no header or terminating end-of-file marker (Hadoop-BAM does not add
        // end-of-file markers), so to merge into a single BAM we concatenate the header with the part files and a
        // terminating end-of-file marker. Note that this has the side effect of being ever-so-slightly less efficient
        // than writing a BAM in one go because the last block of each file isn't completely full.

        Configuration conf = new Configuration();
        FileSystem fs = FileSystem.get(conf);
        String outputParentDir = outputFile.substring(0, outputFile.lastIndexOf("/") + 1);
        // First, check for the _SUCCESS file.
        String successFile = outputFile + "/_SUCCESS";
        if (!fs.exists(new Path(successFile))) {
            throw new GATKException("unable to find " + successFile + " file");
        }
        String tmpPath = outputParentDir + "tmp" + System.currentTimeMillis();
        fs.rename(new Path(outputFile), new Path(tmpPath));
        fs.delete(new Path(outputFile), true);

        try (final OutputStream out = fs.create(new Path(outputFile))) {
            new SAMOutputPreparer().prepareForRecords(out, SAMFormat.BAM, header); // write the header
            mergeInto(out, new Path(tmpPath), conf);
            out.write(BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK); // add the BGZF terminator
        }
    }

    private static void mergeInto(OutputStream out, Path directory, Configuration conf) throws IOException {
        final FileSystem fs = directory.getFileSystem(conf);
        final FileStatus[] parts = fs.globStatus(new Path(directory, "part-r-[0-9][0-9][0-9][0-9][0-9]*"));
        for (final FileStatus part : parts) {
            try (final InputStream in = fs.open(part.getPath())) {
                org.apache.hadoop.io.IOUtils.copyBytes(in, out, conf, false);
            }
        }
        for (final FileStatus part : parts) {
            fs.delete(part.getPath(), false);
        }
    }
}
