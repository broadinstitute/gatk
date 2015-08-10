package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.apache.commons.collections4.iterators.IteratorIterable;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.RecordWriter;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
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
    }

    /**
     * writeReads writes rddReads to outputFile with header as the file header.
     * @param ctx the JavaSparkContext to write.
     * @param outputFile path to the output bam.
     * @param rddReads reads to write.
     * @param header the header to put at the top of the files
     * @param singleFile should the output be a single file or sharded, i.e., outputfile/part-r-00000, outputfile/part-r-00001...
     */
    public static void writeReads(
            final JavaSparkContext ctx, final String outputFile, final JavaRDD<GATKRead> rddReads,
            final SAMFileHeader header, boolean singleFile) throws IOException {
        if (singleFile) {
            writeReadsSingle(ctx, outputFile, rddReads, header);
        } else {
            writeReadsSharded(ctx, outputFile, rddReads, header);

        }

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

        finalOut.saveAsNewAPIHadoopFile(outputFile, GATKRead.class, SAMRecordWritable.class, SparkBAMOutputFormat.class);
        RenameHadoopSingleShard(outputFile);
    }

    private static void RenameHadoopSingleShard(String outputFile) throws IOException {
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