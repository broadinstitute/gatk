package org.broadinstitute.hellbender.engine.spark;

import com.google.common.collect.*;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Lazy;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.spark.Partitioner;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.function.FlatMapFunction2;
import org.apache.spark.api.java.function.MapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.rdd.PartitionCoalescer;
import org.apache.spark.rdd.RDD;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.engine.ShardBoundaryShard;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalLocator;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SerializableSupplier;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import scala.Option;
import scala.Serializable;
import scala.Tuple2;
import scala.reflect.ClassTag;
import scala.reflect.ClassTag$;

import javax.annotation.Nullable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.utils.IntervalUtils.overlaps;

/**
 * Utility methods for sharding {@link Locatable} objects (such as reads) for given intervals, without using a shuffle.
 */
public final class SparkSharder implements AutoCloseable {

    private final JavaSparkContext ctx;
    private final SVIntervalTree<ShardBoundary> shards;
    private final SVIntervalLocator locator;
    private final Lazy<Broadcast<SVIntervalTree<ShardBoundary>>> shardsBroadCast;
    private final Lazy<Broadcast<SVIntervalLocator>> locatorBroadCast;
    private final List<ShardBoundary> shardList;

    /**
     * Creates an {@link SparkSharder} instance.
     * <p>
     *     Subsequent class to instance methods will use the common parameter passed to the constructor:
     *     the {@link JavaSparkContext spark-context}, {@link SAMSequenceDictionary sequence-dictionary},
     *     {@link List<ShardBoundary> shard-intervals} and maximum locatable length.
     * </p>
     *
     * @param ctx the underlying {@link JavaSparkContext}.
     * @param svIntervalLocator Locatable to/from SVInterval translator.
     * @param intervals the intervals to shard.
     */
    public <L extends Locatable> SparkSharder(final JavaSparkContext ctx, final SVIntervalLocator svIntervalLocator,
                        final Collection<L> intervals) {
        this.ctx = Utils.nonNull(ctx);
        this.locator = Utils.nonNull(svIntervalLocator);
        this.shards = new SVIntervalTree<>();
        for (final L element : intervals) {
            this.shards.put(locator.toSVInterval(element), new ShardBoundary(SimpleInterval.valueOf(element)));
        }
        Utils.nonNull(intervals, "the input shard intervals list must be greater than 0");
        this.shardList = Utils.stream(shards).map(SVIntervalTree.Entry::getValue).collect(Collectors.toList());
        this.shardsBroadCast = new Lazy<>(() -> this.ctx.broadcast(shards));
        this.locatorBroadCast = new Lazy<>(() -> this.ctx.broadcast(locator));
    }

    /**
     * Free held up resources.
     */
    @Override
    public void close() {
       if (this.shardsBroadCast.isInitialized()) {
          this.shardsBroadCast.get().destroy(false);
       }
       if (this.locatorBroadCast.isInitialized()) {
          this.locatorBroadCast.get().destroy(false);
       }
    }

   public SparkSharder(final JavaSparkContext ctx, final SAMSequenceDictionary sequenceDictionary,
                       final List<SimpleInterval> intervals, final int shardSize, final int shardPadding) {
       this(ctx, SVIntervalLocator.of(sequenceDictionary), divideIntoShards(sequenceDictionary, intervals, shardSize, shardPadding));
   }

    private static List<ShardBoundary> divideIntoShards(final SAMSequenceDictionary sequenceDictionary, final List<SimpleInterval> intervals, final int shardSize, final int shardPadding) {
        Utils.nonNull(sequenceDictionary);
        ParamUtils.isPositive(shardSize, "the shard size must be positive");
        ParamUtils.isPositiveOrZero(shardPadding, "the shard padding must be 0 or greater");
        final Map<String, List<SimpleInterval>> consolidatedIntervals = IntervalUtils.sortAndMergeIntervals(intervals, true);
        return consolidatedIntervals.values().stream()
                .flatMap(List::stream)
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, shardSize, shardPadding, sequenceDictionary).stream())
                .collect(Collectors.toList());
    }

    /**
     * Creates a RDD of {@link Shard} from an RDD of records.
     * <p>
     *  Each shard contains the records objects that overlap it based on their relevant locations as resolved
     *  by the record to location function provided.
     * </p>
     * <p>
     *     The out {@link ShardRDD} would contain the record at each shard it overlaps exactly once
     *     even if {@code locationOf} would return multiple locations within the same shard.
     * </p>
     *
     * @param elements the RDD of records to shard.
     * @param locationsOf function that provides the relevant locations for each element.
     * @param <L> the {@link Locatable} type.
     * @return an RDD of {@link Shard} of overlapping {@link Locatable} objects (including overlapping only padding)
     */
    public <L> ShardRDD<L> shard(final JavaRDD<L> elements, final SerializableFunction<? super L, List<? extends Locatable>> locationsOf) {
        final Broadcast<SVIntervalTree<ShardBoundary>> shardsBroadCast = this.shardsBroadCast.get();
        final Broadcast<SVIntervalLocator> locatorBroadcast = this.locatorBroadCast.get();
        final JavaPairRDD<ShardBoundary, List<L>> pairRDD = elements.mapPartitionsToPair(it -> {
            final SVIntervalTree<ShardBoundary> shards = shardsBroadCast.getValue();
            final SVIntervalLocator locator = locatorBroadcast.getValue();
            final Map<ShardBoundary, Set<L>> shardsSets = new HashMap<>();
            it.forEachRemaining(l -> {
                for (final Locatable loc : locationsOf.apply(l)) {
                    final Iterator<SVIntervalTree.Entry<ShardBoundary>> overlapping = shards.overlappers(locator.toSVInterval(loc));
                    overlapping.forEachRemaining(entry -> {
                        final Set<L> values = shardsSets.computeIfAbsent(entry.getValue(), o -> new LinkedHashSet<>());
                        values.add(l);
                    });
                }
            });
            return shardsSets.entrySet().stream()
                    .map(entry ->
                            new Tuple2<ShardBoundary, List<L>>(entry.getKey(),
                                    new ArrayList<>(entry.getValue())))
                    .iterator();
        }).reduceByKey((l1, l2) -> {
                    l1.addAll(l2);
                    return l1;
                });
        return new ShardRDD<>(pairRDD);
    }

    /**
     * Creates a RDD of {@link Shard} from an RDD of records.
     * <p>
     *     Each shard contains the records objects that overlap it based on their relevant locations as resolved
     *     by the record to location function provided.
     * </p>
     * <p>
     *     As supposed to {@link #shard(JavaRDD, SerializableFunction) }, this input second order function in this overload allows the composition of the actual location-of function
     *     that reuses computed structures between calls on individual input records.
     * </p>
     *
     * @param elements the RDD that contains the records to shard.
     * @param locationsOfSupplier second order function that provides the function that returns relevant locations for each element.
     * @param <L> the {@link Locatable} type.
     * @return an RDD of {@link Shard} of overlapping {@link Locatable} objects (including overlapping only padding)
     */
    public <L> ShardRDD<L> shard(final JavaRDD<L> elements, final SerializableSupplier<SerializableFunction<? super L, List<? extends Locatable>>> locationsOfSupplier) {
        return shard(elements, Utils.nonNull(locationsOfSupplier).get());
    }

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
     * @return an RDD of {@link Shard} of overlapping {@link Locatable} objects (including overlapping only padding)
     */
    public static <L extends Locatable> JavaRDD<Shard<L>> shard(JavaSparkContext ctx, JavaRDD<L> locatables, Class<L> locatableClass,
                                                                SAMSequenceDictionary sequenceDictionary, List<ShardBoundary> intervals,
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
     * @return an RDD of {@link Shard} of overlapping {@link Locatable} objects (including overlapping only padding)
     */
    public static <L extends Locatable> JavaRDD<Shard<L>> shard(JavaSparkContext ctx, JavaRDD<L> locatables, Class<L> locatableClass,
                                                                SAMSequenceDictionary sequenceDictionary, List<ShardBoundary> intervals,
                                                                int maxLocatableLength, boolean useShuffle) {

        final List<ShardBoundary> paddedIntervals = intervals.stream().map(sb -> new ShardBoundary(sb.getInterval(), sb.getPaddedInterval()) {
            private static final long serialVersionUID = 1L;
            @Override
            public String getContig() {
                return getPaddedInterval().getContig();
            }
            @Override
            public int getStart() {
                return getPaddedInterval().getStart();
            }
            @Override
            public int getEnd() {
                return getPaddedInterval().getEnd();
            }
        }).collect(Collectors.toList());
        if (useShuffle) {
            OverlapDetector<ShardBoundary> overlapDetector = OverlapDetector.create(paddedIntervals);
            Broadcast<OverlapDetector<ShardBoundary>> overlapDetectorBroadcast = ctx.broadcast(overlapDetector);
            JavaPairRDD<ShardBoundary, L> intervalsToLocatables = locatables.flatMapToPair(locatable -> {
                Set<ShardBoundary> overlaps = overlapDetectorBroadcast.getValue().getOverlaps(locatable);
                return overlaps.stream().map(key -> new Tuple2<>(key, locatable)).collect(Collectors.toList()).iterator();
            });
            JavaPairRDD<ShardBoundary, Iterable<L>> grouped = intervalsToLocatables.groupByKey();
            return grouped.map((org.apache.spark.api.java.function.Function<Tuple2<ShardBoundary, Iterable<L>>, Shard<L>>) value -> new ShardBoundaryShard<>(value._1(), value._2()));
        }
        return joinOverlapping(ctx, locatables, locatableClass, sequenceDictionary, paddedIntervals, maxLocatableLength,
                new MapFunction<Tuple2<ShardBoundary, Iterable<L>>, Shard<L>>() {
            private static final long serialVersionUID = 1L;
            @Override
            public Shard<L> call(Tuple2<ShardBoundary, Iterable<L>> value) {
                return new ShardBoundaryShard<>(value._1(), value._2());
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
     * @return never {@code null}.
     */
    private static <L extends Locatable, I extends Locatable, T> JavaRDD<T> joinOverlapping(JavaSparkContext ctx, JavaRDD<L> locatables, Class<L> locatableClass,
                                                                                            SAMSequenceDictionary sequenceDictionary, List<I> intervals,
                                                                                            int maxLocatableLength, MapFunction<Tuple2<I, Iterable<L>>, T> f) {
        return joinOverlapping(ctx, locatables, locatableClass, sequenceDictionary, intervals, maxLocatableLength,
                (FlatMapFunction2<Iterator<L>, Iterator<I>, T>) (locatablesIterator, shardsIterator) -> Iterators.transform(locatablesPerShard(locatablesIterator, shardsIterator, sequenceDictionary, maxLocatableLength), new com.google.common.base.Function<Tuple2<I,Iterable<L>>, T>() {
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
     * @return never {@code null}.
     */
    private static <L extends Locatable, I extends Locatable, T> JavaRDD<T> joinOverlapping(JavaSparkContext ctx, JavaRDD<L> locatables, Class<L> locatableClass,
                                                                                            SAMSequenceDictionary sequenceDictionary, List<I> intervals,
                                                                                            int maxLocatableLength, FlatMapFunction2<Iterator<L>, Iterator<I>, T> f) {

        List<PartitionLocatable<SimpleInterval>> partitionReadExtents = computePartitionReadExtents(locatables, sequenceDictionary, maxLocatableLength);

        // For each interval find which partition it starts and ends in.
        // An interval is processed in the partition it starts in. However, we need to make sure that
        // subsequent partitions are coalesced if needed, so for each partition p find the latest subsequent
        // partition that is needed to read all of the intervals that start in p.
        List<Integer> maxEndPartitionIndexes = new ArrayList<>();
        for (int i = 0; i < locatables.getNumPartitions(); i++) {
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
            indexedIntervals.add(new PartitionLocatable<>(startIndex, interval));
            if (endIndex > maxEndPartitionIndexes.get(startIndex)) {
                maxEndPartitionIndexes.set(startIndex, endIndex);
            }
        }

        JavaRDD<L> coalescedRdd = coalesce(locatables, locatableClass, new RangePartitionCoalescer(maxEndPartitionIndexes));

        // Create an RDD of intervals with the same number of partitions as the locatables, and where each interval
        // is in its start partition.
        JavaRDD<I> intervalsRdd = ctx.parallelize(indexedIntervals)
                .mapToPair(interval ->
                        new Tuple2<>(interval.getPartitionIndex(), interval.getLocatable()))
                .partitionBy(new KeyPartitioner(locatables.getNumPartitions())).values();

        // zipPartitions on coalesced locatable partitions and intervals, and apply the function f
        return coalescedRdd.zipPartitions(intervalsRdd, f);
    }

    /**
     * Turn a pair of iterators over intervals and locatables, into a single iterator over pairs made up of an interval and
     * the locatables that overlap it. Intervals with no overlapping locatables are dropped.
     */
    static <L extends Locatable, I extends Locatable> Iterator<Tuple2<I, Iterable<L>>> locatablesPerShard(
            final Iterator<L> locatables, final Iterator<I> shards, final SAMSequenceDictionary sequenceDictionary,
            final int maxLocatableLength) {
        if (!shards.hasNext()) {
            return Collections.emptyIterator();
        }
        final PeekingIterator<L> peekingLocatables = Iterators.peekingIterator(locatables);
        final PeekingIterator<I> peekingShards = Iterators.peekingIterator(shards);
        final Iterator<Tuple2<I, Iterable<L>>> iterator = new AbstractIterator<Tuple2<I, Iterable<L>>>() {
            // keep track of current and next, since locatables can overlap two shards
            I currentShard = peekingShards.next();
            I nextShard = peekingShards.hasNext() ? peekingShards.next() : null;
            List<L> currentLocatables = Lists.newArrayList();
            List<L> nextLocatables = Lists.newArrayList();

            @Override
            protected Tuple2<I, Iterable<L>> computeNext() {
                if (currentShard == null) {
                    return endOfData();
                }
                while (peekingLocatables.hasNext()) {
                    if (toRightOf(currentShard, peekingLocatables.peek(), sequenceDictionary)) {
                        break;
                    }
                    L locatable = peekingLocatables.next();
                    if (locatable.getContig() != null) {
                        int size = locatable.getEnd() - locatable.getStart() + 1;
                        if (size > maxLocatableLength) {
                            throw new UserException(String.format("Max size of locatable exceeded. Max size is %s, but locatable size is %s. Try increasing shard size and/or padding. Locatable: %s", maxLocatableLength, size, locatable));
                        }
                    }
                    if (overlaps(currentShard, locatable)) {
                        currentLocatables.add(locatable);
                    }
                    if (nextShard != null && overlaps(nextShard, locatable)) {
                        nextLocatables.add(locatable);
                    }
                }
                // current shard is finished, either because the current locatable is to the right of it, or there are no more locatables
                Tuple2<I, Iterable<L>> tuple = new Tuple2<>(currentShard, currentLocatables);
                currentShard = nextShard;
                nextShard = peekingShards.hasNext() ? peekingShards.next() : null;
                currentLocatables = nextLocatables;
                nextLocatables = Lists.newArrayList();
                return tuple;
            }
        };
        return Iterators.filter(iterator, input -> input._2().iterator().hasNext());
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
     * For each partition, find the interval that spans it.
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
                splitPoints.add(new PartitionLocatable<>(i, locatable));
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

    public SerializableFunction<? super Locatable, List<ShardBoundary>> toShardBoundaryFunction() {
        final Broadcast<SVIntervalTree<ShardBoundary>> svIntervalTreeBroadcast = this.shardsBroadCast.get();
        final Broadcast<SVIntervalLocator> locatorBroadcast = this.locatorBroadCast.get();
        return locatable -> {
            final SVIntervalLocator locator = locatorBroadcast.getValue();
            final SVIntervalTree<ShardBoundary> svTree = svIntervalTreeBroadcast.getValue();
            final Iterator<SVIntervalTree.Entry<ShardBoundary>> intervals = svTree.overlappers(locator.toSVInterval(locatable));
            return intervals.hasNext() ? Utils.stream(intervals).map(SVIntervalTree.Entry::getValue).collect(Collectors.toList()) : Collections.emptyList();
        };
    }

    private static class KeyPartitioner extends Partitioner {

        private static final long serialVersionUID = 1L;

        private int numPartitions;

        KeyPartitioner(int numPartitions) {
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

    static class PartitionLocatable<L extends Locatable> implements Locatable, Serializable {
        private static final long serialVersionUID = 1L;

        private final int partitionIndex;
        private final L interval;

        PartitionLocatable(int partitionIndex, L interval) {
            this.partitionIndex = partitionIndex;
            this.interval = interval;
        }

        int getPartitionIndex() {
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
}
