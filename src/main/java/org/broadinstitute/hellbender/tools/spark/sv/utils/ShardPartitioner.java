package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.util.Locatable;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.apache.spark.Partitioner;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.function.Function;

/**
 * Partitioner based on a shard list.
 *
 * Contiguous shards are grouped in partitions trying to keep the same number of shards per partition.
 *
 * <p>
 *     The partition number assigned to a key by the resulting partitioner will "loosely" depend on the shard or shards it overlaps with.
 *     When an object is not overlapped by any shard but maps to a position on a contig with shards present,
 *     its partition would equally be "loosely" dependent of the partitions assigned to shards closest to its position.
 * </p>
 * <p>
 *     When the input key object is not overlapped by any shard, but its contig it amongst the one referred by a shard
 *     in this partitioner, we return a partition number of a shard that is "close" to its location.
 * </p>
 * <p>
 *     When the input object's contig if it is a {@code Locatable} is unknown to the partitioner or it is not a
 *     {@link Locatable} then partitions will be assigned as it was mapping to one of the known contigs
 *     selected at random but consistently based on the contigs name {@link #hashCode}.
 * </p>
 * <p>
 *     Thus this kind of partitioner is guaranteed to do a good job in evenly splitting keys amongst partition if
 *     all keys map to some shard on a known contig and that each shard roughly will contain the same number of
 *     keys. Otherwise, some non-pathological input object key sequences may result in uneven partitioning.
 * </p>
 * <p>
 *     The only guarantee is that for objects A, B, C that map contiguous on the same contig on that order,
 *     if {@code partition(A) == partition(B) == P} then {@code partition(C) == P} also.
 * </p>
 *
 */
public abstract class ShardPartitioner<T> extends Partitioner {

    private final Class<T> clazz;
    private final int numberOfPartitions;
    private final Map<String, Integer> contigToIndex = new HashMap<>();
    private final SVIntervalTree<Integer> partitions = new SVIntervalTree<>();

    /**
     * Creates a partitioner to be used on locatable key objects.
     * @param clazz key object class.
     * @param shards the shards this partitioner will be based on.
     * @param numberOfPartitions number of output partition.
     * @param <T> the class for the key object this partitioner will be used on.
     * @param <L> the type for the shards.
     * @return never {@code null}.
     */
    public static <T extends Locatable, L extends Locatable> ShardPartitioner<T> make(
            final Class<T> clazz, final Collection<L> shards, final int numberOfPartitions) {
        return new ShardPartitioner<T>(clazz, shards, numberOfPartitions) {
            private static final long serialVersionUID = 1L;
            @Override
            public Locatable getLocatable(T key) {
                return key;
            }
        };
    }

    /**
     * Creates a partitioner to be used on non-locatable key objects.
     * @param clazz key object class.
     * @param shards the shards this partitioner will be based on.
     * @param toLocatable function that extracts a {@link Locatable} out of {@link T} typed key objects.
     * @param numberOfPartitions number of output partition.
     * @param <T> the class for the key object this partitioner will be used on.
     * @param <L> the type for the shards.
     * @return never {@code null}.
     */
    public static <T, L extends Locatable> ShardPartitioner<T> make(final Class<T> clazz, final Collection<L> shards,
                                                                    final SerializableFunction<T, Locatable> toLocatable,
                                                                      final int numberOfPartitions) {
        return new ShardPartitioner<T>(clazz, shards, numberOfPartitions) {
            private static final long serialVersionUID = 1L;
            @Override
            public Locatable getLocatable(T key) {
                return toLocatable.apply(key);
            }
        };
    }

    @Override
    public int numPartitions() {
        return numberOfPartitions;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getPartition(final Object key) {
        if (key == null || clazz.isAssignableFrom(key.getClass())) {
            final Locatable locatable = getLocatable(clazz.cast(key));
            return getPartition(locatable);
        } else {
            throw new IllegalArgumentException("key is of the wrong class '" + key.getClass() + "' that does not extends '" + clazz + "'");
        }
    }

    protected abstract Locatable getLocatable(final T key);

    private int getPartition(final Locatable locatable) {
        Utils.nonNull(locatable);
        return getPartition(locatable.getContig(), locatable.getStart(), locatable.getEnd());
    }

    private int getPartition(final String contig, final int start, final int end) {
        final Integer contigIndex = contigToIndex.get(contig);
        final SVInterval query = new SVInterval(contigIndex != null ? contigIndex : Objects.hashCode(contig) % contigToIndex.size(), start, end + 1);
        // By construction minOverlapper won't ever return a null:
        return partitions.minOverlapper(query).getValue();
    }

    /**
     * Create a new shard-partitioner given a sorted collection of shards.
     * @param shards the input shards.
     * @param numberOfPartitions number of partitions for this partitioner.
     * @param <L> the shard type.
     */
    private <L extends Locatable> ShardPartitioner(final Class<T> clazz, final Collection<L> shards, final int numberOfPartitions) {
        Utils.nonNull(clazz);
        Utils.nonNull(shards);
        this.clazz = clazz;
        ParamUtils.isNotEmpty(shards, "shard list");
        ParamUtils.isPositive(numberOfPartitions, "the number of input partitions must be 1 or greater");

        final int numberOfShardsPerPartition = (int) Math.ceil(shards.size() / (float) numberOfPartitions);
        final Iterator<L> shardIterator = shards.iterator();
        String currentContig = null;
        int leftInPartition = numberOfShardsPerPartition; // number of shards left to be included to complete the current
                                                          // partition.
        int currentPartition = 0; // current partition index; first one is 0.
        int currentContigIndex = -1; // current contig index; before the first one is -1.
        int currentPartitionStart = 1; // current partition- start positions set to and open start.
        L previousLocatable = null;
        while (shardIterator.hasNext()) {
            final L locatable = shardIterator.next();
            if (locatable.getContig().equals(currentContig)) {
                if (previousLocatable.getStart() > locatable.getStart() || previousLocatable.getEnd() > locatable.getEnd()) {
                    throw new IllegalArgumentException("the input shard collection contains elements out of order or early elements reach beyond later elements");
                }
                if (leftInPartition == 0) {
                    partitions.put(new SVInterval(currentContigIndex, currentPartitionStart, locatable.getStart()), currentPartition);
                    currentPartitionStart = locatable.getStart();
                }
            } else {
                if (currentContig != null) {
                    partitions.put(new SVInterval(currentContigIndex, currentPartitionStart, Integer.MAX_VALUE), currentPartition);
                }
                currentContigIndex++;
                currentPartitionStart = 1;
                currentContig = locatable.getContig();
                contigToIndex.put(currentContig, currentContigIndex);
            }
            if (leftInPartition-- == 0) {
                leftInPartition = numberOfShardsPerPartition;
                currentPartition++;
            }
            previousLocatable = locatable;
        }
        partitions.put(new SVInterval(currentContigIndex, currentPartitionStart, Integer.MAX_VALUE), currentPartition);
        this.numberOfPartitions = numberOfPartitions;
    }
}
