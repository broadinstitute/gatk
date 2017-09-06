package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.util.Locatable;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.apache.spark.Partitioner;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;

/**
 * Partitioner based a shard list.
 *
 * Contiguous shards are grouped in partitions trying to keep the same number of shards per partition.
 * <p>
 *     This might not always be possible, for example the number of shard in a chromosome is not divisible by
 *     the number of shards in a partition. In these cases, some partitions may contain shards from several chromosomes.
 * </p>
 *
 * <p>
 *     The partition assigned to an object will depend on what shard or shards it overlaps. For {@link Locatable}
 *     instances we use its start position. The object is assigned a partition that could be any of the of the partitions
 *     assigned to any of the shards that overlap that object.
 * </p>
 * <p>
 *     If there are several shards that overlap that object with different partitions assigned to them,
 *     we always consistently return the same partition number for that object.
 * </p>
 * <p>
 *     When the input object is not overlapped by any shard, however its contig it amongst the one referred by a shard
 *     in this partitoner, we return a partition number of a shard that is "close" to its location.
 * </p>
 * <p>
 *     When the input object's contig if it is a {@code Locatable} is unknown to the partitioner or it is not a
 *     {@link Locatable} then a partition will be assigned at random but consistently the same number using the object {@link #hashCode}.
 *     So as long as {@link #hashCode} returns a consistent value, the resulting partion number will be consistently the same
 *     at every invokation.
 * </p>
 */
public final class ShardPartitioner extends Partitioner {

    private static final long serialVersionUID = 1L;

    private final int numberOfPartitions;
    private final Map<String, Integer> contigToIndex = new HashMap<>();
    private final SVIntervalTree<Integer> partitions = new SVIntervalTree<>();

    @Override
    public int numPartitions() {
        return numberOfPartitions;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getPartition(final Object key) {
        if (key instanceof Locatable) {
            return getPartition((Locatable) key);
        } else {
            return key.hashCode() % numberOfPartitions;
        }
    }

    public int getPartition(final Locatable key) {
        Utils.nonNull(key, "unexpected null key");
        return getPartition(key.getContig(), key.getStart());
    }

    public int getPartition(final String contig, final int position) {
        final Integer contigIndex = contigToIndex.get(contig);
        if (contigIndex != null) {
            final SVInterval query = new SVInterval(contigIndex, position - 1, position);
            SVIntervalTree.Entry<Integer> entry = partitions.minOverlapper(query);
            if (entry != null) {
                return entry.getValue();
            }
        }
        return (((((Objects.hashCode(contig) * 47) + position) * 47) + position) * 47) % numberOfPartitions;
    }

    /**
     * Create a new shard-partitioner given a sorted collection of shards.
     * @param shards the input shards.
     * @param numberOfPartitions number of partitions for this partitioner.
     * @param <L> the shard type.
     */
    public <L extends Locatable> ShardPartitioner(final Collection<L> shards, final int numberOfPartitions) {
        Utils.nonNull(shards);
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
                    partitions.put(new SVInterval(currentContigIndex, currentPartitionStart - 1, locatable.getStart() - 1), currentPartition);
                    currentPartitionStart = locatable.getStart();
                }
            } else {
                if (currentContig != null) {
                    partitions.put(new SVInterval(currentContigIndex, currentPartitionStart - 1, Integer.MAX_VALUE), currentPartition);
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
        partitions.put(new SVInterval(currentContigIndex, currentPartitionStart - 1, Integer.MAX_VALUE), currentPartition);
        this.numberOfPartitions = numberOfPartitions;
    }
}
