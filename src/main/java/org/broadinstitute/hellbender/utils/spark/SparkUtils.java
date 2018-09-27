package org.broadinstitute.hellbender.utils.spark;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import htsjdk.samtools.util.RuntimeIOException;
import org.apache.commons.io.FileUtils;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction2;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.*;
import scala.Tuple2;

import java.io.*;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

/**
 * Miscellaneous Spark-related utilities
 */
public final class SparkUtils {
    private static final Logger logger = Logger.getLogger(SparkUtils.class);

    /** Sometimes Spark has trouble destroying a broadcast variable, but we'd like the app to continue anyway. */
    public static <T> void destroyBroadcast(final Broadcast<T> broadcast, final String whatBroadcast ) {
        try {
            broadcast.destroy();
        } catch ( final Exception e ) {
            logger.warn("Failed to destroy broadcast for " + whatBroadcast, e);
        }
    }

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
        final BlockCompressedOutputStream blockCompressedOutputStream = new BlockCompressedOutputStream(outputStream, (File)null);
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
     * Do a total sort of an RDD of {@link GATKRead} according to the sort order in the header.
     * @param reads a JavaRDD of reads which may or may not be sorted
     * @param header a header which specifies the desired new sort order.
     *               Only {@link SAMFileHeader.SortOrder#coordinate} and {@link SAMFileHeader.SortOrder#queryname} are supported.
     *               All others will result in {@link GATKException}
     * @param numReducers number of reducers to use when sorting
     * @return a new JavaRDD or reads which is globally sorted in a way that is consistent with the sort order given in the header
     */
    public static JavaRDD<GATKRead> sortReadsAccordingToHeader(final JavaRDD<GATKRead> reads, final SAMFileHeader header, final int numReducers){
        final SAMFileHeader.SortOrder order = header.getSortOrder();
        switch (order){
            case coordinate:
                return sortUsingElementsAsKeys(reads, new ReadCoordinateComparator(header), numReducers);
            case queryname:
                final JavaRDD<GATKRead> sortedReads = sortUsingElementsAsKeys(reads, new ReadQueryNameComparator(), numReducers);
                return putReadsWithTheSameNameInTheSamePartition(header, sortedReads, JavaSparkContext.fromSparkContext(reads.context()));
            default:
                throw new GATKException("Sort order: " + order + " is not supported.");
        }
    }

    /**
     *   Do a global sort of an RDD using the given comparator.
     *   This method uses the RDD elements themselves as the keys in the spark key/value sort.  This may be inefficient
     *   if the comparator only uses looks at a small fraction of the element to perform the comparison.
     */
    public static <T> JavaRDD<T> sortUsingElementsAsKeys(JavaRDD<T> elements, Comparator<T> comparator, int numReducers) {
        Utils.nonNull(comparator);
        Utils.nonNull(elements);

        // Turn into key-value pairs so we can sort (by key). Values are null so there is no overhead in the amount
        // of data going through the shuffle.
        final JavaPairRDD<T, Void> rddReadPairs = elements.mapToPair(read -> new Tuple2<>(read, (Void) null));

        final JavaPairRDD<T, Void> readVoidPairs;
        if (numReducers > 0) {
            readVoidPairs = rddReadPairs.sortByKey(comparator, true, numReducers);
        } else {
            readVoidPairs = rddReadPairs.sortByKey(comparator);
        }
        return readVoidPairs.keys();
    }

    /**
     * Ensure all reads with the same name appear in the same partition of a queryname sorted RDD.
     * This avoids a global shuffle and only transfers the leading elements from each partition which is fast in most
     * cases.
     *
     * The RDD must be queryname sorted.  If there are so many reads with the same name that they span multiple partitions
     * this will throw {@link GATKException}.
     */
    public static JavaRDD<GATKRead> putReadsWithTheSameNameInTheSamePartition(final SAMFileHeader header, final JavaRDD<GATKRead> reads, final JavaSparkContext ctx) {
        Utils.validateArg(ReadUtils.isReadNameGroupedBam(header), () -> "Reads must be queryname grouped or sorted. " +
                "Actual sort:" + header.getSortOrder() + "  Actual grouping:" +header.getGroupOrder());
        int numPartitions = reads.getNumPartitions();
        final String firstNameInBam = reads.first().getName();
        // Find the first group in each partition
        List<List<GATKRead>> firstReadNameGroupInEachPartition = reads
                .mapPartitions(it -> { PeekingIterator<GATKRead> current = Iterators.peekingIterator(it);
                                List<GATKRead> firstGroup = new ArrayList<>(2);
                                firstGroup.add(current.next());
                                String name = firstGroup.get(0).getName();
                                while (current.hasNext() && current.peek().getName().equals(name)) {
                                    firstGroup.add(current.next());
                                }
                                return Iterators.singletonIterator(firstGroup);
                                })
                .collect();

        // Checking for pathological cases (read name groups that span more than 2 partitions)
        String groupName = null;
        for (List<GATKRead> group : firstReadNameGroupInEachPartition) {
            if (group!=null && !group.isEmpty()) {
                // If a read spans multiple partitions we expect its name to show up multiple times and we don't expect this to work properly
                if (groupName != null && group.get(0).getName().equals(groupName)) {
                    throw new GATKException(String.format("The read name '%s' appeared across multiple partitions this could indicate there was a problem " +
                            "with the sorting or that the rdd has too many partitions, check that the file is queryname sorted and consider decreasing the number of partitions", groupName));
                }
                groupName =  group.get(0).getName();
            }
        }

        // Shift left, so that each partition will be joined with the first read group from the _next_ partition
        List<List<GATKRead>> firstReadInNextPartition = new ArrayList<>(firstReadNameGroupInEachPartition.subList(1, numPartitions));
        firstReadInNextPartition.add(null); // the last partition does not have any reads to add to it

        // Join the reads with the first read from the _next_ partition, then filter out the first reads in this partition
        return reads.zipPartitions(ctx.parallelize(firstReadInNextPartition, numPartitions),
                (FlatMapFunction2<Iterator<GATKRead>, Iterator<List<GATKRead>>, GATKRead>) (it1, it2) -> {
            PeekingIterator<GATKRead> current = Iterators.peekingIterator(it1);
            String firstName = current.peek().getName();
            // Make sure we don't remove reads from the first partition
            if (!firstNameInBam.equals(firstName)) {
                // skip the first read name group in the _current_ partition since it will be handled in the previous partition
                while (current.hasNext() && current.peek() != null && current.peek().getName().equals(firstName)) {
                    current.next();
                }
            }
            // append the first reads in the _next_ partition to the _current_ partition
            PeekingIterator<List<GATKRead>> next = Iterators.peekingIterator(it2);
            if (next.hasNext() && next.peek() != null) {
                return Iterators.concat(current, next.peek().iterator());
            }
            return current;
        });
    }
}
