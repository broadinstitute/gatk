package org.broadinstitute.hellbender.utils.spark;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.*;
import org.apache.commons.io.FileUtils;
import org.apache.spark.Partitioner;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.function.FlatMapFunction2;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.rdd.RDD;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import scala.Tuple2;
import scala.reflect.ClassTag;
import scala.reflect.ClassTag$;


import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

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
     * Join an RDD of reads with a set of intervals, and apply a function to process the reads that overlap each interval.
     * @param ctx the Spark Context
     * @param reads the reads RDD, must be coordinate sorted
     * @param readClass the class of the reads, must be a subclass of {@link Locatable}
     * @param intervals the collection of intervals to apply the function to
     * @param f the function to process intervals and overlapping reads with
     * @param <R> the read type
     * @param <I> the interval type
     * @param <T> the return type of <code>f</code>
     * @return
     */
    public static <R extends Locatable, I extends Locatable, T> JavaRDD<T> joinOverlapping(JavaSparkContext ctx, JavaRDD<R> reads, Class<R> readClass,
                                                                                           List<I> intervals,
                                                                                           FlatMapFunction2<Iterator<R>, Iterator<I>, T> f) {
        // Find the total extent of all reads in each partition. This requires that the input reads
        // are scanned over (to find end points). This should be faster than shuffling though.
        List<PartitionLocatable<SimpleInterval>> partitionReadExtents = reads.mapPartitionsWithIndex(new Function2<Integer, Iterator<R>, Iterator<PartitionLocatable<SimpleInterval>>>() {
            private static final long serialVersionUID = 1L;

            @Override
            public Iterator<PartitionLocatable<SimpleInterval>> call(Integer index, Iterator<R> it) throws Exception {
                if (!it.hasNext()) {
                    return Collections.emptyIterator();
                }
                List<PartitionLocatable<SimpleInterval>> extents = new ArrayList<>();
                R read = it.next();
                String contig = read.getContig();
                int minStart = read.getStart();
                int maxEnd = read.getEnd();
                while (it.hasNext()) {
                    R next = it.next();
                    if (!contig.equals(next.getContig())) {
                        extents.add(new PartitionLocatable<>(index, new SimpleInterval(contig, minStart, maxEnd)));
                        contig = next.getContig();
                        minStart = next.getStart();
                        maxEnd = next.getEnd();
                        continue;
                    }
                    maxEnd = Math.max(maxEnd, next.getEnd());
                }
                extents.add(new PartitionLocatable<>(index, new SimpleInterval(contig, minStart, maxEnd)));
                return extents.iterator();
            }
        }, true).collect();

        // For each interval find which partition it starts and ends in.
        // An interval is processed in the partition it starts in. However, we need to make sure that
        // subsequent partitions are coalesced if needed, so for each partition p find the latest subsequent
        // partition that is needed to read all of the intervals that start in p.
        List<Integer> maxEndPartitionIndexes = new ArrayList<>();
        for (int i = 0; i < reads.getNumPartitions(); i++) {
            maxEndPartitionIndexes.add(i);
        }
        OverlapDetector<PartitionLocatable<SimpleInterval>> overlapDetector = OverlapDetector.create(partitionReadExtents);
        List<PartitionLocatable<I>> indexedIntervals = new ArrayList<>();
        for (I interval : intervals) {
            int[] partitionIndexes = overlapDetector.getOverlaps(interval).stream()
                    .mapToInt(PartitionLocatable::getPartitionIndex).toArray();
            if (partitionIndexes.length == 0) {
                // interval does not overlap any partition - skip it
                continue;
            }
            Arrays.sort(partitionIndexes);
            int startIndex = partitionIndexes[0];
            int endIndex = partitionIndexes[partitionIndexes.length - 1];
            indexedIntervals.add(new PartitionLocatable<I>(startIndex, interval));
            if (endIndex > maxEndPartitionIndexes.get(startIndex)) {
                maxEndPartitionIndexes.set(startIndex, endIndex);
            }
        }

        JavaRDD<R> coalescedRdd = coalesce(reads, readClass, new RangePartitionCoalescer(maxEndPartitionIndexes));

        // Create an RDD of intervals with the same number of partitions as the reads, and where each interval
        // is in its start partition.
        JavaRDD<I> intervalsRdd = ctx.parallelize(indexedIntervals)
                .mapToPair(interval ->
                        new Tuple2<>(interval.getPartitionIndex(), interval.getLocatable()))
                .partitionBy(new KeyPartitioner(reads.getNumPartitions())).values();

        // zipPartitions on coalesced read partitions and intervals, and apply the function f
        return coalescedRdd.zipPartitions(intervalsRdd, f);
    }

    private static <T> JavaRDD<T> coalesce(JavaRDD<T> rdd, Class<T> cls, PartitionCoalescer partitionCoalescer) {
        RDD<T> coalescedRdd = new CoalescedRDD<>(rdd.rdd(), rdd.getNumPartitions(), partitionCoalescer, cls);
        ClassTag<T> tag = ClassTag$.MODULE$.apply(cls);
        return new JavaRDD<>(coalescedRdd, tag);
    }

    private static class KeyPartitioner extends Partitioner {

        private static final long serialVersionUID = 1L;

        private int numPartitions;

        public KeyPartitioner(int numPartitions) {
            this.numPartitions = numPartitions;
        }

        @Override
        public int numPartitions() {
            return numPartitions;
        }

        @Override
        public int getPartition(Object key) {
            return (Integer) key;
        }

    }

    private static class PartitionLocatable<L extends Locatable> implements Locatable {
        private static final long serialVersionUID = 1L;

        private final int partitionIndex;
        private final L interval;


        public PartitionLocatable(int partitionIndex, L interval) {
            this.partitionIndex = partitionIndex;
            this.interval = interval;
        }

        public int getPartitionIndex() {
            return partitionIndex;
        }

        public L getLocatable() {
            return interval;
        }

        @Override
        public String getContig() {
            return interval.getContig();
        }

        @Override
        public int getStart() {
            return interval.getStart();
        }

        @Override
        public int getEnd() {
            return interval.getEnd();
        }

        @Override
        public String toString() {
            return "PartitionLocatable{" +
                    "partitionIndex=" + partitionIndex +
                    ", interval='" + interval + '\'' +
                    '}';
        }
    }

    public static <R extends Locatable, I extends Locatable, T> JavaRDD<T> joinOverlappingNaive(JavaSparkContext ctx, JavaRDD<R> reads, Class<R> readClass,
                                                                                           List<I> intervals,
                                                                                                FlatMapFunction<Tuple2<I,Iterable<R>>, T> f) {
        OverlapDetector<I> overlapDetector = OverlapDetector.create(intervals);
        Broadcast<OverlapDetector<I>> overlapDetectorBroadcast = ctx.broadcast(overlapDetector);
        JavaPairRDD<I, R> intervalsToReads = reads.flatMapToPair(read -> {
            Set<I> overlaps = overlapDetectorBroadcast.getValue().getOverlaps(read);
            return overlaps.stream().map(key -> new Tuple2<>(key, read)).collect(Collectors.toList());
        });
        JavaPairRDD<I, Iterable<R>> grouped = intervalsToReads.groupByKey();
        return grouped.flatMap(f);
    }
}
