package org.broadinstitute.hellbender.utils.spark;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import htsjdk.samtools.util.RuntimeIOException;
import org.apache.commons.io.FileUtils;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.fs.FileSystem;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;


import java.io.*;
import java.util.Comparator;

/**
 * Miscellaneous Spark-related utilities
 */
public final class SparkUtils {

    private SparkUtils() {}

    /**
     * Converts a headerless Hadoop bam shard (eg., a part0000, part0001, etc. file produced by
     * {@link org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink}) into a readable bam file
     * by adding a header and a BGZF terminator.
     *
     * This method is not intended for use with Hadoop bam shards that already have a header -- these shards are
     * already readable using samtools. Currently {@link ReadsSparkSink} saves the "shards" with a header for the
     * {@link ReadsWriteFormat#SHARDED} case, and without a header for the {@link ReadsWriteFormat#SINGLE} case.
     *
     * @param bamShard The headerless Hadoop bam shard to convert
     * @param header header for the BAM file to be created
     * @param destination path to which to write the new BAM file
     */
    public static void convertHeaderlessHadoopBamShardToBam( final File bamShard, final SAMFileHeader header, final File destination ) {
        try ( FileOutputStream outStream = new FileOutputStream(destination) ) {
            writeBAMHeaderToStream(header, outStream);
            FileUtils.copyFile(bamShard, outStream);
            outStream.write(BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK);
        }
        catch ( IOException e ) {
            throw new UserException("Error writing to " + destination.getAbsolutePath(), e);
        }
    }

    /**
     * Private helper method for {@link #convertHeaderlessHadoopBamShardToBam} that takes a SAMFileHeader and writes it
     * to the provided `OutputStream`, correctly encoded for the BAM format and preceded by the BAM magic bytes.
     *
     * @param samFileHeader SAM header to write
     * @param outputStream stream to write the SAM header to
     */
    private static void writeBAMHeaderToStream( final SAMFileHeader samFileHeader, final OutputStream outputStream ) {
        final BlockCompressedOutputStream blockCompressedOutputStream = new BlockCompressedOutputStream(outputStream, null);
        final BinaryCodec outputBinaryCodec = new BinaryCodec(new DataOutputStream(blockCompressedOutputStream));

        final String headerString;
        final Writer stringWriter = new StringWriter();
        new SAMTextHeaderCodec().encode(stringWriter, samFileHeader, true);
        headerString = stringWriter.toString();

        outputBinaryCodec.writeBytes(ReadUtils.BAM_MAGIC);

        // calculate and write the length of the SAM file header text and the header text
        outputBinaryCodec.writeString(headerString, true, false);

        // write the sequences binarily.  This is redundant with the text header
        outputBinaryCodec.writeInt(samFileHeader.getSequenceDictionary().size());
        for (final SAMSequenceRecord sequenceRecord: samFileHeader.getSequenceDictionary().getSequences()) {
            outputBinaryCodec.writeString(sequenceRecord.getSequenceName(), true, true);
            outputBinaryCodec.writeInt(sequenceRecord.getSequenceLength());
        }

        try {
            blockCompressedOutputStream.flush();
        } catch (final IOException ioe) {
            throw new RuntimeIOException(ioe);
        }
    }

    /**
     * Determine if the <code>targetPath</code> exists.
     * @param ctx JavaSparkContext
     * @param targetPath the <code>org.apache.hadoop.fs.Path</code> object to check
     * @return true if the targetPath exists, otherwise false
     */
    public static boolean pathExists(final JavaSparkContext ctx, final Path targetPath) {
        Utils.nonNull(ctx);
        Utils.nonNull(targetPath);
        try {
            final FileSystem fs = targetPath.getFileSystem(ctx.hadoopConfiguration());
            return fs.exists(targetPath);
        } catch (IOException e) {
            throw new UserException("Error validating existence of path " + targetPath + ": " + e.getMessage());
        }
    }

    /**
     * Sorts the given reads in coordinate sort order.
     * @param reads the reads to sort
     * @param header the reads header, which must specify coordinate sort order
     * @param numReducers the number of reducers to use; a value of 0 means use the default number of reducers
     * @return a sorted RDD of reads
     */
    public static JavaRDD<GATKRead> coordinateSortReads(final JavaRDD<GATKRead> reads, final SAMFileHeader header, final int numReducers) {
        Utils.validate(header.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate), "Header must specify coordinate sort order, but was" + header.getSortOrder());

        // Turn into key-value pairs so we can sort (by key). Values are null so there is no overhead in the amount
        // of data going through the shuffle.
        final JavaPairRDD<GATKRead, Void> rddReadPairs = reads.mapToPair(read -> new Tuple2<>(read, (Void) null));

        // do a total sort so that all the reads in partition i are less than those in partition i+1
        final Comparator<GATKRead> comparator = new ReadCoordinateComparator(header);
        final JavaPairRDD<GATKRead, Void> readVoidPairs;
        if (numReducers > 0) {
            readVoidPairs = rddReadPairs.sortByKey(comparator, true, numReducers);
        } else {
            readVoidPairs = rddReadPairs.sortByKey(comparator);
        }
        return readVoidPairs.keys();
    }

    /**
     * Sorts the given reads according to the sort order in the header.
     * @param reads the reads to sort
     * @param header the header specifying the sort order
     * @param numReducers the number of reducers to use; a vlue of 0 means use the default number of reducers
     * @return a sorted RDD of reads
     */
    public static JavaRDD<SAMRecord> sortReads(final JavaRDD<SAMRecord> reads, final SAMFileHeader header, final int numReducers) {
        // Turn into key-value pairs so we can sort (by key). Values are null so there is no overhead in the amount
        // of data going through the shuffle.
        final JavaPairRDD<SAMRecord, Void> rddReadPairs = reads.mapToPair(read -> new Tuple2<>(read, (Void) null));

        // do a total sort so that all the reads in partition i are less than those in partition i+1
        final Comparator<SAMRecord> comparator = getSAMRecordComparator(header);
        final JavaPairRDD<SAMRecord, Void> readVoidPairs;
        if (comparator == null){
            readVoidPairs = rddReadPairs; //no sort
        } else if (numReducers > 0) {
            readVoidPairs = rddReadPairs.sortByKey(comparator, true, numReducers);
        } else {
            readVoidPairs = rddReadPairs.sortByKey(comparator);
        }
        return readVoidPairs.keys();
    }

    //Returns the comparator to use or null if no sorting is required.
    private static Comparator<SAMRecord> getSAMRecordComparator(final SAMFileHeader header) {
        switch (header.getSortOrder()){
            case coordinate: return new HeaderlessSAMRecordCoordinateComparator(header);
            case duplicate:
            case queryname:
            case unsorted:   return header.getSortOrder().getComparatorInstance();
            default:         return null; //NOTE: javac warns if you have this (useless) default BUT it errors out if you remove this default.
        }
    }
}
