package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.apache.commons.collections4.iterators.IteratorIterable;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapred.FileAlreadyExistsException;
import org.apache.hadoop.mapreduce.JobContext;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.RecordWriter;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.parquet.avro.AvroParquetOutputFormat;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.bdgenomics.adam.models.RecordGroupDictionary;
import org.bdgenomics.adam.models.SequenceDictionary;
import org.bdgenomics.formats.avro.AlignmentRecord;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GATKReadToBDGAlignmentRecordConverter;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.seqdoop.hadoop_bam.KeyIgnoringBAMOutputFormat;
import org.seqdoop.hadoop_bam.SAMRecordWritable;
import scala.Tuple2;

import java.io.IOException;

/**
 * ReadsSparkSink writes GATKReads to a file. This code lifts from the HadoopGenomics/Hadoop-BAM
 * read writing code as well as from bigdatagenomics/adam.
 */
public class ReadsSparkSink {
    // TODO: Make ReadsSparkSink able to also write sharded BAMs (#854).

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

    /**
     * writeReads writes rddReads to outputFile with header as the file header.
     * @param ctx the JavaSparkContext to write.
     * @param outputFile path to the output bam.
     * @param rddReads reads to write.
     * @param header the header to put at the top of the files
     * @param format should the output be a single file, sharded, ADAM, etc.
     */
    public static void writeReads(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<GATKRead> rddReads,
            final SAMFileHeader header, ReadsWriteFormat format) throws IOException {
        if (format.equals(ReadsWriteFormat.SINGLE)) {
            writeReadsSingle(ctx, outputFile, rddReads, header);
        } else if (format.equals(ReadsWriteFormat.SHARDED)) {
            writeReadsSharded(ctx, outputFile, rddReads, header);
        } else if (format.equals(ReadsWriteFormat.ADAM)) {
            writeReadsADAM(ctx, outputFile, rddReads, header);
        }
    }

    private static void writeReadsADAM(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<GATKRead> rddGATKReads,
            final SAMFileHeader header) throws IOException {
        SequenceDictionary seqDict = SequenceDictionary.fromSAMSequenceDictionary(header.getSequenceDictionary());
        RecordGroupDictionary readGroups = RecordGroupDictionary.fromSAMHeader(header);
        JavaPairRDD<Void, AlignmentRecord> rddAlignmentRecords =
                rddGATKReads.map(read -> GATKReadToBDGAlignmentRecordConverter.convert(read, header, seqDict, readGroups))
                        .mapToPair(alignmentRecord -> new Tuple2<>(null, alignmentRecord));
        // instantiating a Job is necessary here in order to set the Hadoop Configuration...
        Job job = Job.getInstance(ctx.hadoopConfiguration());
        // ...here, which sets a config property that the AvroParquetOutputFormat needs when writing data. Specifically,
        // we are writing the Avro schema to the Configuration as a JSON string. The AvroParquetOutputFormat class knows
        // how to translate objects in the Avro data model to the Parquet primitives that get written.
        AvroParquetOutputFormat.setSchema(job, AlignmentRecord.getClassSchema());
        rddAlignmentRecords.saveAsNewAPIHadoopFile(
                outputFile, Void.class, AlignmentRecord.class, AvroParquetOutputFormat.class, job.getConfiguration());
    }

    private static void writeReadsSharded(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<GATKRead> rddReads,
            final SAMFileHeader header) {
        // Set the header on the main thread.
        SparkBAMOutputFormat.setHeader(header);
        // MyOutputFormat is a static class, so we need to copy the header to each worker then call
        final Broadcast<SAMFileHeader> broadcast = ctx.broadcast(header);
        JavaRDD<GATKRead> gatkReadJavaRDD = rddReads.mapPartitions(gatkReadIterator -> {
            SparkBAMOutputFormat.setHeader(broadcast.getValue());
            return new IteratorIterable<>(gatkReadIterator);
        });

        // The expected format for writing is JavaPairRDD where the key is ignored and the value is a
        // SAMRecordWritable. We use GATKRead as the key, so we can sort the reads before writing them to a file.
        JavaPairRDD<GATKRead, SAMRecordWritable> rddSamRecordWriteable = gatkReadJavaRDD.mapToPair(gatkRead -> {
            SAMRecord samRecord = gatkRead.convertToSAMRecord(broadcast.getValue());
            SAMRecordWritable samRecordWritable = new SAMRecordWritable();
            samRecordWritable.set(samRecord);
            return new Tuple2<>(gatkRead, samRecordWritable);
        });

        rddSamRecordWriteable.saveAsNewAPIHadoopFile(outputFile, GATKRead.class, SAMRecordWritable.class, SparkBAMOutputFormat.class);
    }

    private static void writeReadsSingle(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<GATKRead> rddReads,
            final SAMFileHeader header) throws IOException {
        // Set the header on the main thread.
        SparkBAMOutputFormat.setHeader(header);

        // The expected format for writing is JavaPairRDD where the key is ignored and the value is a
        // SAMRecordWritable. We use GATKRead as the key, so we can sort the reads before writing them to a file.
        JavaPairRDD<GATKRead, SAMRecordWritable> rddSamRecordWriteable = rddReads.mapToPair(gatkRead -> {
            SAMRecord samRecord = gatkRead.convertToSAMRecord(header);
            SAMRecordWritable samRecordWritable = new SAMRecordWritable();
            samRecordWritable.set(samRecord);
            return new Tuple2<>(gatkRead, samRecordWritable);
        });

        // coalesce is needed for writing to a single bam.
        final JavaPairRDD<GATKRead, SAMRecordWritable> out =
                rddSamRecordWriteable.sortByKey(new ReadCoordinateComparator(header)).coalesce(1);

        // MyOutputFormat is a static class, so we need to copy the header to each worker then call
        // MyOutputFormat.setHeader.
        final Broadcast<SAMFileHeader> broadcast = ctx.broadcast(header);

        final JavaPairRDD<GATKRead, SAMRecordWritable> finalOut = out.mapPartitions(tuple2Iterator -> {
            SparkBAMOutputFormat.setHeader(broadcast.getValue());
            return new IteratorIterable<>(tuple2Iterator);
        }).mapToPair(t -> t);

        deleteHadoopFile(outputFile);
        finalOut.saveAsNewAPIHadoopFile(outputFile, GATKRead.class, SAMRecordWritable.class, SparkBAMOutputFormat.class);
        renameHadoopSingleShard(outputFile);
    }

    private static void deleteHadoopFile(String fileToObliterate) throws IOException {
        Configuration conf = new Configuration();
        FileSystem fs = FileSystem.get(conf);
        fs.delete(new Path(fileToObliterate),true);
    }

    private static void renameHadoopSingleShard(String outputFile) throws IOException {
        // At this point, the file is in outputFile/part-r-0000. We move it out of that dir and into it's own well-named
        // file.  This will not work if the last step was a merge (assuming saveAsNewAPIHadoopFile doesn't cause the
        // last step to be a reduce).
        Configuration conf = new Configuration();
        FileSystem fs = FileSystem.get(conf);
        String ouputParentDir = outputFile.substring(0, outputFile.lastIndexOf("/") + 1);
        // First, check for the _SUCCESS file.
        String successFile = outputFile + "/_SUCCESS";
        if (!fs.exists(new Path(successFile))) {
            throw new GATKException("unable to find " + successFile + " file");
        }
        String tmpPath = ouputParentDir + "tmp" + System.currentTimeMillis();
        fs.rename(new Path(outputFile + "/part-r-00000"), new Path(tmpPath));
        fs.delete(new Path(outputFile), true);
        fs.rename(new Path(tmpPath), new Path(outputFile));
    }
}
