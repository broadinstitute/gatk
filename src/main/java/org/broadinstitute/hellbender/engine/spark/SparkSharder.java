package org.broadinstitute.hellbender.engine.spark;

import com.google.common.base.Function;
import com.google.common.collect.*;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.spark.Partitioner;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.*;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.rdd.PartitionCoalescer;
import org.apache.spark.rdd.RDD;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Option;
import scala.Tuple2;
import scala.reflect.ClassTag;
import scala.reflect.ClassTag$;

import javax.annotation.Nullable;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.utils.IntervalUtils.overlaps;

/**
 * Utility methods for sharding {@link Locatable} objects (such as reads) for given intervals, without using a shuffle.
 */
public class SparkSharder {
    /**
     * Create an RDD of {@link Shard} from an RDD of coordinate sorted {@link Locatable} <i>without using a shuffle</i>.
     * Each shard contains the {@link Locatable} objects that overlap it (including overlapping only padding).
     * @param ctx the Spark Context
     * @param locatables the RDD of {@link Locatable}, must be coordinate sorted
     * @param locatableClass the class of the {@link Locatable} objects in the RDD
     * @param sequenceDictionary the sequence dictionary to use to find contig lengths
     * @param intervals the {@link ShardBoundary} objects to create shards for, must be coordinate sorted
     * @param maxLocatableLength the maximum length of a {@link Locatable}, if any is larger than this size then an exception will be thrown
     * @param <L> the {@link Locatable} type
     * @param <SB> the {@link ShardBoundary} type
     * @return an RDD of {@link Shard} of overlapping {@link Locatable} objects (including overlapping only padding)
     */
    public static <L extends Locatable, SB extends ShardBoundary> JavaRDD<Shard<L>> shard(JavaSparkContext ctx, JavaRDD<L> locatables, Class<L> locatableClass,
                                                                SAMSequenceDictionary sequenceDictionary, List<SB> intervals,
                                                                int maxLocatableLength) {
        return shard(ctx, locatables, locatableClass, sequenceDictionary, intervals, maxLocatableLength, false);
    }

    /**
     * Create an RDD of {@link Shard} from an RDD of coordinate sorted {@link Locatable} <i>without using a shuffle</i>,
     * and where the intervals for shards are specified as an RDD, rather than a list.
     * Each shard contains the {@link Locatable} objects that overlap it (including overlapping only padding).
     * @param ctx the Spark Context
     * @param locatables the RDD of {@link Locatable}, must be coordinate sorted
     * @param locatableClass the class of the {@link Locatable} objects in the RDD
     * @param sequenceDictionary the sequence dictionary to use to find contig lengths
     * @param intervals the {@link ShardBoundary} objects to create shards for, must be coordinate sorted
     * @param maxLocatableLength the maximum length of a {@link Locatable}, if any is larger than this size then an exception will be thrown
     * @param <L> the {@link Locatable} type
     * @param <SB> the {@link ShardBoundary} type
     * @return an RDD of {@link Shard} of overlapping {@link Locatable} objects (including overlapping only padding)
     */
    public static <L extends Locatable, SB extends ShardBoundary> JavaRDD<Shard<L>> shard(JavaSparkContext ctx, JavaRDD<L> locatables, Class<L> locatableClass,
                                                                SAMSequenceDictionary sequenceDictionary, JavaRDD<SB> intervals,
                                                                int maxLocatableLength) {
        return shard(ctx, locatables, locatableClass, sequenceDictionary, intervals, maxLocatableLength, false);
    }

    /**
     * Create an RDD of {@link Shard} from an RDD of coordinate sorted {@link Locatable}, optionally using a shuffle.
     * A shuffle is typically only needed for correctness testing, since it usually has a significant performance impact.
     * @param ctx the Spark Context
     * @param locatables the RDD of {@link Locatable}, must be coordinate sorted
     * @param locatableClass the class of the {@link Locatable} objects in the RDD
     * @param sequenceDictionary the sequence dictionary to use to find contig lengths
     * @param intervals the {@link ShardBoundary} objects to create shards for, must be coordinate sorted
     * @param maxLocatableLength the maximum length of a {@link Locatable}, if any is larger than this size then an exception will be thrown
     * @param useShuffle whether to use a shuffle or not
     * @param <L> the {@link Locatable} type
     * @param <SB> the {@link ShardBoundary} type
     * @return an RDD of {@link Shard} of overlapping {@link Locatable} objects (including overlapping only padding)
     */
    public static <L extends Locatable, SB extends ShardBoundary> JavaRDD<Shard<L>> shard(JavaSparkContext ctx, JavaRDD<L> locatables, Class<L> locatableClass,
                                                                SAMSequenceDictionary sequenceDictionary, List<SB> intervals,
                                                                int maxLocatableLength, boolean useShuffle) {

        List<ShardBoundary> paddedIntervals = intervals.stream().map(ShardBoundary::paddedShardBoundary).collect(Collectors.toList());
        if (useShuffle) {
            OverlapDetector<ShardBoundary> overlapDetector = OverlapDetector.create(paddedIntervals);
            Broadcast<OverlapDetector<ShardBoundary>> overlapDetectorBroadcast = ctx.broadcast(overlapDetector);
            JavaPairRDD<ShardBoundary, L> intervalsToLocatables = locatables.flatMapToPair(locatable -> {
                Set<ShardBoundary> overlaps = overlapDetectorBroadcast.getValue().getOverlaps(locatable);
                return overlaps.stream().map(key -> new Tuple2<>(key, locatable)).collect(Collectors.toList()).iterator();
            });
            JavaPairRDD<ShardBoundary, Iterable<L>> grouped = intervalsToLocatables.groupByKey();
            return grouped.map((org.apache.spark.api.java.function.Function<Tuple2<ShardBoundary, Iterable<L>>, Shard<L>>) value -> value._1().createShard(value._2()));
        }
        return joinOverlapping(ctx, locatables, locatableClass, sequenceDictionary, paddedIntervals, maxLocatableLength,
                new MapFunction<Tuple2<ShardBoundary, Iterable<L>>, Shard<L>>() {
            private static final long serialVersionUID = 1L;
            @Override
            public Shard<L> call(Tuple2<ShardBoundary, Iterable<L>> value) {
                return value._1().createShard(value._2());
            }
        });
    }

    private static <L extends Locatable, SB extends ShardBoundary> JavaRDD<Shard<L>> shard(JavaSparkContext ctx, JavaRDD<L> locatables, Class<L> locatableClass,
                                                                SAMSequenceDictionary sequenceDictionary, JavaRDD<SB> intervals,
                                                                int maxLocatableLength, boolean useShuffle) {

        JavaRDD<ShardBoundary> paddedIntervals = intervals.map(ShardBoundary::paddedShardBoundary);
        if (useShuffle) {
            throw new UnsupportedOperationException("Shuffle not supported when sharding an RDD of intervals.");
        }
        return joinOverlapping(ctx, locatables, locatableClass, sequenceDictionary, paddedIntervals, maxLocatableLength,
                new MapFunction<Tuple2<ShardBoundary, Iterable<L>>, Shard<L>>() {
                    private static final long serialVersionUID = 1L;
                    @Override
                    public Shard<L> call(Tuple2<ShardBoundary, Iterable<L>> value) {
                        return value._1().createShard(value._2());
                    }
                });
    }

    /**
     * Join an RDD of locatables with a set of intervals, and apply a function to process the locatables that overlap each interval.
     * @param ctx the Spark Context
     * @param locatables the locatables RDD, must be coordinate sorted
     * @param locatableClass the class of the locatables, must be a subclass of {@link Locatable}
     * @param sequenceDictionary the sequence dictionary to use to find contig lengths
     * @param intervals the collection of intervals to apply the function to
     * @param maxLocatableLength the maximum length of a {@link Locatable}, if any is larger than this size then an exception will be thrown
     * @param f the function to process intervals and overlapping locatables with
     * @param <L> the {@link Locatable} type
     * @param <I> the interval type
     * @param <T> the return type of <code>f</code>
     * @return
     */
    private static <L extends Locatable, I extends Locatable, T> JavaRDD<T> joinOverlapping(JavaSparkContext ctx, JavaRDD<L> locatables, Class<L> locatableClass,
                                                                                            SAMSequenceDictionary sequenceDictionary, List<I> intervals,
                                                                                            int maxLocatableLength, MapFunction<Tuple2<I, Iterable<L>>, T> f) {
        return joinOverlapping(ctx, locatables, locatableClass, sequenceDictionary, intervals, maxLocatableLength,
                (FlatMapFunction2<Iterator<L>, Iterator<I>, T>) (locatablesIterator, shardsIterator) -> Iterators.transform(locatablesPerShard(locatablesIterator, shardsIterator, sequenceDictionary, maxLocatableLength), new Function<Tuple2<I,Iterable<L>>, T>() {
                    @Nullable
                    @Override
                    public T apply(@Nullable Tuple2<I, Iterable<L>> input) {
                        try {
                            return f.call(input);
                        } catch (Exception e) {
                            throw new RuntimeException(e);
                        }
                    }
                }));
    }

    private static <L extends Locatable, I extends Locatable, T> JavaRDD<T> joinOverlapping(JavaSparkContext ctx, JavaRDD<L> locatables, Class<L> locatableClass,
                                                                                            SAMSequenceDictionary sequenceDictionary, JavaRDD<I> intervals,
                                                                                            int maxLocatableLength, MapFunction<Tuple2<I, Iterable<L>>, T> f) {
        return joinOverlapping(ctx, locatables, locatableClass, sequenceDictionary, intervals, maxLocatableLength,
                (FlatMapFunction2<Iterator<L>, Iterator<I>, T>) (locatablesIterator, shardsIterator) -> Iterators.transform(locatablesPerShard(locatablesIterator, shardsIterator, sequenceDictionary, maxLocatableLength), new Function<Tuple2<I,Iterable<L>>, T>() {
                    @Nullable
                    @Override
                    public T apply(@Nullable Tuple2<I, Iterable<L>> input) {
                        try {
                            return f.call(input);
                        } catch (Exception e) {
                            throw new RuntimeException(e);
                        }
                    }
                }));
    }

    /**
     * Join an RDD of locatables with a set of intervals, and apply a function to process the locatables that overlap each interval.
     * This differs from {@link #joinOverlapping(JavaSparkContext, JavaRDD, Class, SAMSequenceDictionary, List, int, MapFunction)}
     * in that the function to apply is given two iterators: one over intervals, and one over locatables (for the partition),
     * and it is up to the function implemention to find overlaps between intervals and locatables.
     * @param ctx the Spark Context
     * @param locatables the locatables RDD, must be coordinate sorted
     * @param locatableClass the class of the locatables, must be a subclass of {@link Locatable}
     * @param sequenceDictionary the sequence dictionary to use to find contig lengths
     * @param intervals the collection of intervals to apply the function to
     * @param maxLocatableLength the maximum length of a {@link Locatable}, if any is larger than this size then an exception will be thrown
     * @param f the function to process intervals and overlapping locatables with
     * @param <L> the {@link Locatable} type
     * @param <I> the interval type
     * @param <T> the return type of <code>f</code>
     * @return
     */
    private static <L extends Locatable, I extends Locatable, T> JavaRDD<T> joinOverlapping(JavaSparkContext ctx, JavaRDD<L> locatables, Class<L> locatableClass,
                                                                                            SAMSequenceDictionary sequenceDictionary, List<I> intervals,
                                                                                            int maxLocatableLength, FlatMapFunction2<Iterator<L>, Iterator<I>, T> f) {

        return joinOverlapping(ctx, locatables, locatableClass, sequenceDictionary, ctx.parallelize(intervals), maxLocatableLength, f);
    }

    private static <L extends Locatable, I extends Locatable, T> JavaRDD<T> joinOverlapping(JavaSparkContext ctx, JavaRDD<L> locatables, Class<L> locatableClass,
                                                                                            SAMSequenceDictionary sequenceDictionary, JavaRDD<I> intervals,
                                                                                            int maxLocatableLength, FlatMapFunction2<Iterator<L>, Iterator<I>, T> f) {

        List<PartitionLocatable<SimpleInterval>> partitionReadExtents = computePartitionReadExtents(locatables, sequenceDictionary, maxLocatableLength);
        List<SimpleInterval> firstLocatablesList = partitionReadExtents.stream().map(PartitionLocatable::getLocatable).collect(Collectors.toList());
        Broadcast<List<SimpleInterval>> firstLocatablesBroadcast = ctx.broadcast(firstLocatablesList);

        // For each interval find which partition it starts and ends in.
        // An interval is processed in the partition it starts in. However, we need to make sure that
        // subsequent partitions are coalesced if needed, so for each partition p find the latest subsequent
        // partition that is needed to read all of the intervals that start in p.
        OverlapDetector<PartitionLocatable<SimpleInterval>> overlapDetector = OverlapDetector.create(partitionReadExtents);
        Broadcast<OverlapDetector<PartitionLocatable<SimpleInterval>>> overlapDetectorBroadcast = ctx.broadcast(overlapDetector);
        JavaRDD<PartitionLocatable<I>> indexedIntervals = intervals.map(interval -> {
            int[] partitionIndexes = overlapDetectorBroadcast.getValue().getOverlaps(interval).stream()
                    .mapToInt(PartitionLocatable::getPartitionIndex).toArray();
            if (partitionIndexes.length == 0) {
                final List<SimpleInterval> firstLocatables = firstLocatablesBroadcast.getValue();
                // interval does not overlap any partition - add it to the one after the interval start
                int i = Collections.binarySearch(firstLocatables, new SimpleInterval(interval), (o1, o2) -> IntervalUtils.compareLocatables(o1, o2, sequenceDictionary));
                if (i >= 0) {
                    throw new IllegalStateException(); // TODO: no overlaps, yet start of interval matches a partition read extent start
                }
                int insertionPoint = -i - 1;
                if (insertionPoint == firstLocatables.size()) {
                    insertionPoint = firstLocatables.size() - 1;
                }
                return new PartitionLocatable<>(insertionPoint, interval);
            }
            Arrays.sort(partitionIndexes);
            int startIndex = partitionIndexes[0];
            int endIndex = partitionIndexes[partitionIndexes.length - 1];
            return new PartitionLocatable<>(startIndex, endIndex, interval);
        });

        // Create an RDD of intervals with the same number of partitions as the locatables, and where each interval
        // is in its start partition. Within each partition, intervals are sorted by IntervalUtils#compareLocatables.
        JavaRDD<PartitionLocatable<I>> indexedIntervalsRepartitioned = indexedIntervals
                .mapToPair(interval ->
                        new Tuple2<>(interval, (Void) null))
                .repartitionAndSortWithinPartitions(new PartitionLocatablePartitioner(locatables.getNumPartitions()), new PartitionLocatableComparator<I>(sequenceDictionary))
                .keys();

        indexedIntervalsRepartitioned.cache(); // cache since we need to do two calculations on the intervals

        // Find the end partition index for each partition.
        Map<Integer, Integer> maxEndPartitionIndexesMap = indexedIntervalsRepartitioned.mapToPair((PairFunction<PartitionLocatable<I>, Integer, Integer>) partitionLocatable ->
                new Tuple2<>(partitionLocatable.getPartitionIndex(), partitionLocatable.getEndPartitionIndex()))
                .reduceByKey((Function2<Integer, Integer, Integer>) Math::max)
                .collectAsMap();
        List<Integer> maxEndPartitionIndexes = IntStream.range(0, locatables.getNumPartitions()).boxed().collect(Collectors.toList());
        maxEndPartitionIndexesMap.forEach((startIndex, endIndex) -> {
            if (endIndex > maxEndPartitionIndexes.get(startIndex)) {
                maxEndPartitionIndexes.set(startIndex, endIndex);
            }
        });

        JavaRDD<L> coalescedRdd = coalesce(locatables, locatableClass, new RangePartitionCoalescer(maxEndPartitionIndexes));

        // zipPartitions on coalesced locatable partitions and intervals, and apply the function f
        return coalescedRdd.zipPartitions(indexedIntervalsRepartitioned.map(PartitionLocatable::getLocatable), f);
    }

    /**
     * Turn a pair of iterators over intervals and locatables, into a single iterator over pairs made up of an interval and
     * the locatables that overlap it. Intervals with no overlapping locatables are included.
     */
    static <L extends Locatable, I extends Locatable> Iterator<Tuple2<I, Iterable<L>>> locatablesPerShard(Iterator<L> locatables, Iterator<I> shards, SAMSequenceDictionary sequenceDictionary, int maxLocatableLength) {
        if (!shards.hasNext()) {
            return Collections.emptyIterator();
        }
        PeekingIterator<I> peekingShards = Iterators.peekingIterator(shards);
        Iterator<Tuple2<I, Iterable<L>>> iterator = new AbstractIterator<Tuple2<I, Iterable<L>>>() {
            Queue<PendingShard<L, I>> pendingShards = new ArrayDeque<>();

            @Override
            protected Tuple2<I, Iterable<L>> computeNext() {
                Tuple2<I, Iterable<L>> nextShard = null;

                while (locatables.hasNext() && nextShard == null) {
                    L locatable = locatables.next();
                    if (locatable.getContig() != null) {
                        int size = locatable.getEnd() - locatable.getStart() + 1;
                        if (size > maxLocatableLength) {
                            throw new UserException(String.format("Max size of locatable exceeded. Max size is %s, but locatable size is %s. Try increasing shard size and/or padding. Locatable: %s", maxLocatableLength, size, locatable));
                        }
                    }

                    // Add any new shards that start before the end of the read to the queue
                    while (peekingShards.hasNext() && !IntervalUtils.isAfter(peekingShards.peek(), locatable, sequenceDictionary)) {
                        pendingShards.add(new PendingShard<>(peekingShards.next()));
                    }

                    // Add the read to any shards that it overlaps
                    for (PendingShard<L, I> pendingShard : pendingShards) {
                        if (overlaps(pendingShard, locatable)) {
                            pendingShard.addLocatable(locatable);
                        }
                    }

                    // A pending shard only becomes ready once our reads iterator has advanced beyond the end of its extended span
                    // (this ensures that we've loaded all reads that belong in the new shard)
                    if (!pendingShards.isEmpty() && IntervalUtils.isAfter(locatable, pendingShards.peek(), sequenceDictionary)) {
                        nextShard = pendingShards.poll().get();
                    }
                }

                if (!locatables.hasNext()) {
                    // Pull on intervals until it is exhausted
                    while (peekingShards.hasNext()) {
                        pendingShards.add(new PendingShard<>(peekingShards.next()));
                    }

                    // Grab the next pending shard if there is one, unless we already have a shard ready to go
                    if (!pendingShards.isEmpty() && nextShard == null) {
                        nextShard = pendingShards.poll().get();
                    }
                }

                if (nextShard == null) {
                    return endOfData();
                }

                return nextShard;
            }
        };
        return iterator;
    }

    private static class PendingShard<L extends Locatable, I extends Locatable> implements Locatable {
        private I interval;
        private List<L> locatables = new ArrayList<>();

        public PendingShard(I interval) {
            this.interval = interval;
        }

        public void addLocatable(L locatable) {
            locatables.add(locatable);
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

        public Tuple2<I, Iterable<L>> get() {
            return new Tuple2<>(interval, locatables);
        }
    }

    /**
     * @return <code>true</code> if the locatable is to the right of the given interval
     */
    private static <I extends Locatable, L extends Locatable> boolean toRightOf(I interval, L locatable, SAMSequenceDictionary sequenceDictionary) {
        int intervalContigIndex = sequenceDictionary.getSequenceIndex(interval.getContig());
        int locatableContigIndex = sequenceDictionary.getSequenceIndex(locatable.getContig());
        return (intervalContigIndex == locatableContigIndex && interval.getEnd() < locatable.getStart()) // locatable on same contig, to the right
                || intervalContigIndex < locatableContigIndex; // locatable on subsequent contig
    }

    /**
     * For each partition, find the interval that spans it, ordered by start position.
     */
    static <L extends Locatable> List<PartitionLocatable<SimpleInterval>> computePartitionReadExtents(JavaRDD<L> locatables, SAMSequenceDictionary sequenceDictionary, int maxLocatableLength) {
        // Find the first locatable in each partition. This is very efficient since only the first record in each partition is read.
        // If a partition is empty then set the locatable to null
        List<PartitionLocatable<L>> allSplitPoints = locatables.mapPartitions(
                (FlatMapFunction<Iterator<L>, PartitionLocatable<L>>) it -> ImmutableList.of(new PartitionLocatable<>(-1, it.hasNext() ? it.next() : null)).iterator()
        ).collect();
        List<PartitionLocatable<L>> splitPoints = new ArrayList<>(); // fill in index and remove nulls (empty partitions)
        for (int i = 0; i < allSplitPoints.size(); i++) {
            L locatable = allSplitPoints.get(i).getLocatable();
            if (locatable != null) {
                splitPoints.add(new PartitionLocatable<L>(i, locatable));
            }
        }
        List<PartitionLocatable<SimpleInterval>> extents = new ArrayList<>();
        for (int i = 0; i < splitPoints.size(); i++) {
            PartitionLocatable<L> splitPoint = splitPoints.get(i);
            int partitionIndex = splitPoint.getPartitionIndex();
            Locatable current = splitPoint.getLocatable();
            int intervalContigIndex = sequenceDictionary.getSequenceIndex(current.getContig());
            Utils.validate(intervalContigIndex != -1, "Contig not found in sequence dictionary: " + current.getContig());
            final Locatable next;
            final int nextContigIndex;
            if (i < splitPoints.size() - 1) {
                next = splitPoints.get(i + 1);
                nextContigIndex = sequenceDictionary.getSequenceIndex(next.getContig());
                Utils.validate(nextContigIndex != -1, "Contig not found in sequence dictionary: " + next.getContig());
            } else {
                next = null;
                nextContigIndex = sequenceDictionary.getSequences().size();
            }
            if (intervalContigIndex == nextContigIndex) { // same contig
                addPartitionReadExtent(extents, partitionIndex, current.getContig(), current.getStart(), next.getStart() + maxLocatableLength);
            } else {
                // complete current contig
                SAMSequenceRecord seq = sequenceDictionary.getSequence(current.getContig());
                Utils.validate(seq != null, "Contig not found in sequence dictionary: " + current.getContig());
                int contigEnd = seq.getSequenceLength();
                addPartitionReadExtent(extents, partitionIndex, current.getContig(), current.getStart(), contigEnd);
                // add any whole contigs up to next (exclusive)
                for (int contigIndex = intervalContigIndex + 1; contigIndex < nextContigIndex; contigIndex++) {
                    SAMSequenceRecord sequence = sequenceDictionary.getSequence(contigIndex);
                    Utils.validate(sequence != null, "Contig index not found in sequence dictionary: " + contigIndex);
                    addPartitionReadExtent(extents, partitionIndex, sequence.getSequenceName(), 1, sequence.getSequenceLength());
                }
                // add start of next contig
                if (next != null) {
                    addPartitionReadExtent(extents, partitionIndex, next.getContig(), 1, next.getStart() + maxLocatableLength);
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
        RDD<T> coalescedRdd = rdd.rdd().coalesce(rdd.getNumPartitions(), false, Option.apply(partitionCoalescer), null);
        ClassTag<T> tag = ClassTag$.MODULE$.apply(cls);
        return new JavaRDD<>(coalescedRdd, tag);
    }

    /**
     * Assigns {@link PartitionLocatable} objects to their start partition.
     */
    private static class PartitionLocatablePartitioner extends Partitioner {
        private static final long serialVersionUID = 1L;

        private int numPartitions;

        public PartitionLocatablePartitioner(int numPartitions) {
            this.numPartitions = numPartitions;
        }

        @Override
        public int numPartitions() {
            return numPartitions;
        }

        @Override
        public int getPartition(Object key) {
            return ((PartitionLocatable) key).getPartitionIndex();
        }
    }

    /**
     * Compares {@link PartitionLocatable} objects using a {@link htsjdk.samtools.SAMSequenceDictionary} sequence ordering.
     * @param <L> the interval type
     */
    private static class PartitionLocatableComparator<L extends Locatable> implements Comparator<PartitionLocatable<L>>, Serializable {
        private static final long serialVersionUID = 1L;
        private final SAMSequenceDictionary sequenceDictionary;

        private PartitionLocatableComparator(SAMSequenceDictionary sequenceDictionary) {
            this.sequenceDictionary = sequenceDictionary;
        }

        @Override
        public int compare(PartitionLocatable<L> pl1, PartitionLocatable<L> pl2) {
            return IntervalUtils.compareLocatables(pl1.getLocatable(), pl2.getLocatable(), this.sequenceDictionary);
        }
    }

    /**
     * Encapsulates the start and end partitions for an interval.
     * @param <L> the interval type
     */
    static class PartitionLocatable<L extends Locatable> implements Locatable {
        private static final long serialVersionUID = 1L;

        private final int partitionIndex;
        private final int endPartitionIndex;
        private final L interval;

        public PartitionLocatable(int partitionIndex, L interval) {
            this(partitionIndex, partitionIndex, interval);
        }

        public PartitionLocatable(int partitionIndex, int endPartitionIndex, L interval) {
            this.partitionIndex = partitionIndex;
            this.endPartitionIndex = endPartitionIndex;
            this.interval = interval;
        }

        public int getPartitionIndex() {
            return partitionIndex;
        }

        public int getEndPartitionIndex() {
            return endPartitionIndex;
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
}
