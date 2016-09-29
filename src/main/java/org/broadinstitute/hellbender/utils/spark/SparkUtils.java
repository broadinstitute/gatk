package org.broadinstitute.hellbender.utils.spark;

import com.google.common.base.Function;
import com.google.common.base.Preconditions;
import com.google.common.collect.*;
import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterables;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.*;
import org.apache.commons.io.FileUtils;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.fs.FileSystem;
import org.apache.spark.Partitioner;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.function.FlatMapFunction2;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.rdd.RDD;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceSource;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import scala.Tuple2;
import scala.reflect.ClassTag;
import scala.reflect.ClassTag$;


import javax.annotation.Nullable;
import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.utils.IntervalUtils.overlaps;

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
                                                                                           SAMSequenceDictionary sequenceDictionary, List<I> intervals,
                                                                                           FlatMapFunction2<Iterator<R>, Iterator<I>, T> f) {

        List<PartitionLocatable<SimpleInterval>> partitionReadExtents = computePartitionReadExtents(reads, sequenceDictionary);

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

    public static <R extends Locatable, I extends Locatable, T> JavaRDD<T> joinOverlapping(JavaSparkContext ctx, JavaRDD<R> reads, Class<R> readClass,
                                                                                           SAMSequenceDictionary sequenceDictionary, List<I> intervals,
                                                                                           FlatMapFunction<Tuple2<I,Iterable<R>>, T> f) {
        return joinOverlapping(ctx, reads, readClass, sequenceDictionary, intervals,
                (FlatMapFunction2<Iterator<R>, Iterator<I>, T>) (iterator, iIterator) -> {
                    Iterator<Iterator<T>> transform = Iterators.transform(readsPerShard(iterator, iIterator), new Function<Tuple2<I, Iterable<R>>, Iterator<T>>() {
                        @Nullable
                        @Override
                        public Iterator<T> apply(@Nullable Tuple2<I, Iterable<R>> input) {
                            try {
                                return f.call(input).iterator();
                            } catch (Exception e) {
                                throw new RuntimeException(e);
                            }
                        }
                    });
                    return () -> Iterators.concat(transform);
                });
    }

    public static <T> JavaPairRDD<GATKRead, Tuple2<T, ReferenceBases>> addBases(JavaSparkContext ctx, final ReferenceMultiSource referenceSource,
                                                                                final JavaPairRDD<GATKRead, T> keyedByRead,
                                                                                final SAMSequenceDictionary sequenceDictionary) {
        List<PartitionLocatable<SimpleInterval>> partitionReadExtents = computePartitionReadExtents(keyedByRead.map(t -> t._1()), sequenceDictionary);

        Broadcast<ReferenceMultiSource> bReferenceSource = ctx.broadcast(referenceSource);
        JavaRDD<ReferenceBases> referenceBasesRDD = ctx.parallelize(partitionReadExtents)
                .mapToPair(interval ->
                        new Tuple2<>(interval.getPartitionIndex(), interval.getLocatable()))
                .partitionBy(new KeyPartitioner(keyedByRead.getNumPartitions()))
                .values()
                .map(interval -> bReferenceSource.getValue().getReferenceBases(null, interval)); // get bases for partition on executor

        return keyedByRead.zipPartitions(referenceBasesRDD, (FlatMapFunction2<Iterator<Tuple2<GATKRead, T>>, Iterator<ReferenceBases>, Tuple2<GATKRead, Tuple2<T, ReferenceBases>>>) (readTupleIterator, referenceBasesIterator) -> () -> {
            ReferenceBases refBases = Iterators.getOnlyElement(referenceBasesIterator);
            return Iterators.transform(readTupleIterator, new Function<Tuple2<GATKRead, T>, Tuple2<GATKRead, Tuple2<T, ReferenceBases>>>() {
                @Nullable
                @Override
                public Tuple2<GATKRead, Tuple2<T, ReferenceBases>> apply(@Nullable Tuple2<GATKRead, T> input) {
                    return new Tuple2<>(input._1(), new Tuple2<>(input._2(), refBases));
                }
            });
        }).mapToPair(t -> t);
    }

    public static <R extends Locatable, I extends Locatable> Iterator<Tuple2<I, Iterable<R>>> readsPerShard(Iterator<R> reads, Iterator<I> shards) {
        PeekingIterator<R> peekingReads = Iterators.peekingIterator(reads);
        PeekingIterator<I> peekingShards = Iterators.peekingIterator(shards);
        return new AbstractIterator<Tuple2<I, Iterable<R>>>() {
            // keep track of current and next, since reads can overlap two shards
            I currentShard = peekingShards.next();
            I nextShard = peekingShards.hasNext() ? peekingShards.next() : null;
            List<R> currentReads = Lists.newArrayList();
            List<R> nextReads = Lists.newArrayList();
            @Override
            protected Tuple2<I, Iterable<R>> computeNext() {
                if (currentShard == null) {
                    return endOfData();
                }
                while (peekingReads.hasNext()) {
                    if (toRightOf(currentShard, peekingReads.peek())) {
                        break;
                    }
                    R read = peekingReads.next();
                    if (overlaps(currentShard, read)) {
                        currentReads.add(read);
                    }
                    if (nextShard != null && overlaps(nextShard, read)) {
                        nextReads.add(read);
                    }
                }
                // current shard is finished, either because the current read is to the right of it, or there are no more reads
                Tuple2<I, Iterable<R>> tuple = new Tuple2<>(currentShard, currentReads);
                currentShard = nextShard;
                nextShard = peekingShards.hasNext() ? peekingShards.next() : null;
                currentReads = nextReads;
                nextReads = Lists.newArrayList();
                return tuple;
            }
        };
    }

    private static <I extends Locatable, R extends Locatable> boolean toRightOf(I currentShard, R read) {
        return currentShard.getEnd() < read.getStart();
    }

    static <R extends Locatable> List<PartitionLocatable<SimpleInterval>> computePartitionReadExtents(JavaRDD<R> reads, SAMSequenceDictionary sequenceDictionary) {
        // Find the first read in each partition. This is very efficient since only the first record in each partition is read.
        List<R> splitPoints = reads.mapPartitions((FlatMapFunction<Iterator<R>, R>) it -> ImmutableList.of(it.next())).collect();
        List<PartitionLocatable<SimpleInterval>> extents = new ArrayList<>();
        for (int i = 0; i < splitPoints.size(); i++) {
            Locatable current = splitPoints.get(i);
            int intervalContigIndex = sequenceDictionary.getSequenceIndex(current.getContig());
            final Locatable next;
            final int nextContigIndex;
            if (i < splitPoints.size() - 1) {
                next = splitPoints.get(i + 1);
                nextContigIndex = sequenceDictionary.getSequenceIndex(next.getContig());
            } else {
                next = null;
                nextContigIndex = sequenceDictionary.getSequences().size();
            }
            if (intervalContigIndex == nextContigIndex) { // same contig
                addPartitionReadExtent(extents, i, current.getContig(), current.getStart(), next.getEnd()); // assumes reads are same size
            } else {
                // complete current contig
                int contigEnd = sequenceDictionary.getSequence(current.getContig()).getSequenceLength();
                addPartitionReadExtent(extents, i, current.getContig(), current.getStart(), contigEnd);
                // add any whole contigs up to next (exclusive)
                for (int contigIndex = intervalContigIndex + 1; contigIndex < nextContigIndex; contigIndex++) {
                    SAMSequenceRecord sequence = sequenceDictionary.getSequence(contigIndex);
                    addPartitionReadExtent(extents, i, sequence.getSequenceName(), 1, sequence.getSequenceLength());
                }
                // add start of next contig
                if (next != null) {
                    addPartitionReadExtent(extents, i, next.getContig(), 1, next.getEnd()); // assumes reads are same size
                }
            }
        }
        return extents;
    }

    private static void addPartitionReadExtent(List<PartitionLocatable<SimpleInterval>> extents, int partitionIndex, String contig, int start, int end) {
        SimpleInterval extent = new SimpleInterval(contig, start, end);
        extents.add(new PartitionLocatable<>(partitionIndex, extent));
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

    static class PartitionLocatable<L extends Locatable> implements Locatable {
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

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            PartitionLocatable<?> that = (PartitionLocatable<?>) o;

            if (partitionIndex != that.partitionIndex) return false;
            return interval.equals(that.interval);

        }

        @Override
        public int hashCode() {
            int result = partitionIndex;
            result = 31 * result + interval.hashCode();
            return result;
        }
    }

    public static <R extends Locatable, I extends Locatable, T> JavaRDD<T> joinOverlappingShuffle(JavaSparkContext ctx, JavaRDD<R> reads, Class<R> readClass,
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
