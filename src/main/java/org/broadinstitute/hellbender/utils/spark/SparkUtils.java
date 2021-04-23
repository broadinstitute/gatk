package org.broadinstitute.hellbender.utils.spark;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
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
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.*;
import scala.Tuple2;

import java.io.*;
import java.net.URI;
import java.util.*;

/**
 * Miscellaneous Spark-related utilities
 */
public final class SparkUtils {
    private static final Logger logger = LogManager.getLogger(SparkUtils.class);

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
     * @param targetURI the <code>org.apache.hadoop.fs.Path</code> URI to check
     * @return true if the targetPath exists, otherwise false
     */
    public static boolean hadoopPathExists(final JavaSparkContext ctx, final URI targetURI) {
        Utils.nonNull(ctx);
        Utils.nonNull(targetURI);
        try {
            final Path targetHadoopPath = new Path(targetURI);
            final FileSystem fs = targetHadoopPath.getFileSystem(ctx.hadoopConfiguration());
            return fs.exists(targetHadoopPath);
        } catch (IOException e) {
            throw new UserException("Error validating existence of path " + targetURI + ": " + e.getMessage());
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
    public static JavaRDD<GATKRead> putReadsWithTheSameNameInTheSamePartition( final SAMFileHeader header,
                                                                               final JavaRDD<GATKRead> reads,
                                                                               final JavaSparkContext ctx ) {
        Utils.validateArg(ReadUtils.isReadNameGroupedBam(header), () -> "Reads must be queryname grouped or sorted. " +
                "Actual sort:" + header.getSortOrder() + "  Actual grouping:" + header.getGroupOrder());

        // Find the first group in each partition
        final List<List<GATKRead>> firstReadNameGroupInEachPartition = reads
                .mapPartitions(it -> {
                    if ( !it.hasNext() ) {
                        return Iterators.singletonIterator(Collections.<GATKRead>emptyList());
                    }
                    final List<GATKRead> firstGroup = new ArrayList<>(2);
                    final GATKRead firstRead = it.next();
                    firstGroup.add(firstRead);
                    final String groupName = firstRead.getName();
                    while ( it.hasNext() ) {
                        final GATKRead read = it.next();
                        if ( !groupName.equals(read.getName()) ) {
                            break;
                        }
                        firstGroup.add(read);
                    }
                    return Iterators.singletonIterator(firstGroup);
                })
                .collect();

        // Shift left, so that each partition will be zipped with the first read group from the _next_ partition
        final int numPartitions = reads.getNumPartitions();
        final List<List<GATKRead>> firstGroupFromNextPartition =
                new ArrayList<>(firstReadNameGroupInEachPartition.subList(1, numPartitions));
        firstGroupFromNextPartition.add(Collections.emptyList()); // the last partition does not have any reads to add to it

        // Take care of the situation where an entire partition contains reads with the same name
        // (unlikely, but could happen with very long reads, or very small partitions).
        for ( int idx = numPartitions - 1; idx >= 1; --idx ) {
            final List<GATKRead> curGroup = firstGroupFromNextPartition.get(idx);
            if ( !curGroup.isEmpty() ) {
                final String groupName = curGroup.get(0).getName();
                int idx2 = idx;
                while ( --idx2 >= 0 ) {
                    final List<GATKRead> prevGroup = firstGroupFromNextPartition.get(idx2);
                    if ( !prevGroup.isEmpty() ) {
                        if ( groupName.equals(prevGroup.get(0).getName()) ) {
                            prevGroup.addAll(curGroup);
                            curGroup.clear();
                        }
                        break;
                    }
                }
            }
        }

        // Peel off the first group in each partition
        final int[] firstGroupSizes = firstReadNameGroupInEachPartition.stream().mapToInt(List::size).toArray();
        firstGroupSizes[0] = 0; // first partition has no predecessor to handle its first group of reads
        JavaRDD<GATKRead> readsSansFirstGroup = reads.mapPartitionsWithIndex( (idx, itr) ->
            { int groupSize = firstGroupSizes[idx];
              while ( itr.hasNext() && groupSize-- > 0 ) {
                  itr.next();
              }
              return itr; }, true);

        // Zip up the remaining reads with the first read group from the _next_ partition
        return readsSansFirstGroup.zipPartitions(ctx.parallelize(firstGroupFromNextPartition, numPartitions),
                (it1, it2) -> Iterators.concat(it1, it2.next().iterator()));
    }

    /**
     * Like <code>groupByKey</code>, but assumes that values are already sorted by key, so no shuffle is needed,
     * which is much faster.
     * @param rdd the input RDD
     * @param <K> type of keys
     * @param <V> type of values
     * @return an RDD where each the values for each key are grouped into an iterable collection
     */
    public static <K, V> JavaPairRDD<K, Iterable<V>> spanByKey(JavaPairRDD<K, V> rdd) {
        return rdd.mapPartitionsToPair(SparkUtils::getSpanningIterator);
    }

    /**
     * An iterator that groups values having the same key into iterable collections.
     * @param iterator an iterator over key-value pairs
     * @param <K> type of keys
     * @param <V> type of values
     * @return an iterator over pairs of keys and grouped values
     */
    public static <K, V> Iterator<Tuple2<K, Iterable<V>>> getSpanningIterator(Iterator<Tuple2<K, V>> iterator) {
        final PeekingIterator<Tuple2<K, V>> iter = Iterators.peekingIterator(iterator);
        return new AbstractIterator<Tuple2<K, Iterable<V>>>() {
            @Override
            protected Tuple2<K, Iterable<V>> computeNext() {
                K key = null;
                List<V> group = Lists.newArrayList();
                while (iter.hasNext()) {
                    if (key == null) {
                        Tuple2<K, V> next = iter.next();
                        key = next._1();
                        V value = next._2();
                        group.add(value);
                        continue;
                    }
                    K nextKey = iter.peek()._1(); // don't advance...
                    if (nextKey.equals(key)) {
                        group.add(iter.next()._2()); // .. unless the keys match
                    } else {
                        return new Tuple2<>(key, group);
                    }
                }
                if (key != null) {
                    return new Tuple2<>(key, group);
                }
                return endOfData();
            }
        };
    }

    /**
     * Sort reads into queryname order if they are not already sorted
     */
    public static JavaRDD<GATKRead> querynameSortReadsIfNecessary(JavaRDD<GATKRead> reads, int numReducers, SAMFileHeader header) {
        JavaRDD<GATKRead> sortedReadsForMarking;
        if (ReadUtils.isReadNameGroupedBam(header)) {
            sortedReadsForMarking = reads;
        } else {
            header.setSortOrder(SAMFileHeader.SortOrder.queryname);
            sortedReadsForMarking = sortReadsAccordingToHeader(reads, header, numReducers);
        }
        return sortedReadsForMarking;
    }
}
