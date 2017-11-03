package org.broadinstitute.hellbender.tools.spark.sv.utils;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.util.Locatable;
import org.apache.spark.Partitioner;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.Serializable;
import java.util.*;

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
//@DefaultSerializer(ShardPartitioner.Serializer.class)
public final class ShardPartitioner<T extends Locatable> extends Partitioner {

    private static final long serialVersionUID = 1L;

    private final Class<T> clazz;
    private final int numberOfPartitions;
    private final Map<String, Integer> contigToIndex = new LinkedHashMap<>();
    private final SVIntervalTree<Integer> partitions;

    @SuppressWarnings("unchecked")
    private ShardPartitioner(final Serialized serialized) {
        try {
            clazz = (Class<T>) Class.forName(serialized.clazzName);
        } catch (final ClassNotFoundException ex) {
            throw new GATKException.ShouldNeverReachHereException(ex);
        }
        numberOfPartitions = serialized.numberOfPartitions;
        for (int i = 0; i < serialized.contigNames.length; i++) {
            contigToIndex.put(serialized.contigNames[i], serialized.contigIndex[i]);
        }
        partitions = new SVIntervalTree<>();
        for (int i = 0; i < serialized.intervalStart.length; i++) {
            partitions.put(new SVInterval(serialized.intervalIndex[i], serialized.intervalStart[i], serialized.intervalEnd[i]),
                    serialized.intervalPartition[i]);
        }
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
            final Locatable locatable = clazz.cast(key);
            return getPartition(locatable);
        } else {
            throw new IllegalArgumentException("key is of the wrong class '" + key.getClass() + "' that does not extends '" + clazz + "'");
        }
    }

    private int getPartition(final Locatable locatable) {
        Utils.nonNull(locatable);
        return getPartition(locatable.getContig(), locatable.getStart(), locatable.getEnd());
    }

    private int getPartition(final String contig, final int start, final int end) {
        final Integer contigIndex = contigToIndex.get(contig);
        final SVInterval query = new SVInterval(contigIndex != null ? contigIndex : Objects.hashCode(contig) % contigToIndex.size(), start, end + 1);
        // By construction minOverlapper won't ever return a null:
        //return 0;
        final SVIntervalTree.Entry<Integer> minOverlapper = partitions.minOverlapper(query);
        return minOverlapper.getValue();
    }

    @SuppressWarnings("unchecked")
    private ShardPartitioner(final Kryo kryo, final Input input) {
        final String clazzName = input.readString();
        try {
            this.clazz = (Class<T>) Class.forName(clazzName);
        } catch (final ClassNotFoundException ex) {
            throw new GATKException("unknown Locatable class " + clazzName);
        }
        if (!Locatable.class.isAssignableFrom(clazz)) {
            throw new GATKException("the element class does not extend Locatable: " + clazzName);
        }
        this.numberOfPartitions = input.readInt();
        final int numberOfContigs = contigToIndex.size();
        for (int i = 0; i < numberOfContigs; i++) {
            contigToIndex.put(input.readString(), input.readInt());
        }
        partitions = (SVIntervalTree<Integer>) kryo.getSerializer(SVIntervalTree.class).read(kryo, input, (Class<SVIntervalTree<Integer>>) (Class) SVIntervalTree.class);
    }

    /**
     * Create a new shard-partitioner given a sorted collection of shards.
     * @param shards the input shards.
     * @param numberOfPartitions number of partitions for this partitioner.
     * @param <L> the shard type.
     */
    public <L extends Locatable> ShardPartitioner(final Class<T> clazz, final Collection<L> shards, final int numberOfPartitions) {
        Utils.nonNull(clazz);
        Utils.nonNull(shards);
        this.clazz = clazz;
        ParamUtils.isNotEmpty(shards, "shard list");
        ParamUtils.isPositive(numberOfPartitions, "the number of input partitions must be 1 or greater");
        final int numberOfShardsPerPartition = (int) Math.ceil(shards.size() / (float) numberOfPartitions);
        final Iterator<L> shardIterator = shards.iterator();
        partitions = new SVIntervalTree<>();

        String currentContig = null;
        int leftInPartition = numberOfShardsPerPartition; // number of shards left to be included to complete the current
                                                          // partition.
        int currentPartition = 0; // current partition index; first one is 0.
        int currentContigIndex = -1; // current contig index; before the first one is -1.
        int currentPartitionStart = 0; // current partition- start positions set to and open start.
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

    @SuppressWarnings("unchecked")
    private void serialize(final Kryo kryo, final Output output) {
        output.writeString(clazz.getName());
        output.writeInt(numberOfPartitions);
        output.writeInt(contigToIndex.size());
        for (final Map.Entry<String, Integer> entry : contigToIndex.entrySet()) {
            output.writeString(entry.getKey());
            output.writeInt(entry.getValue());
        }
        kryo.getSerializer(partitions.getClass()).write(kryo, output, partitions);
    }

    public static final class Serializer<T extends Locatable> extends com.esotericsoftware.kryo.Serializer<ShardPartitioner<T>> {
        @Override
        public void write(final Kryo kryo, final Output output, final ShardPartitioner<T> interval ) {
            interval.serialize(kryo, output);
        }

        @Override
        public  ShardPartitioner<T> read(final Kryo kryo, final Input input, final Class<ShardPartitioner<T>> klass ) {
            return new ShardPartitioner<>(kryo, input);
        }
    }

    // Needed for serialization since {@link #partitions SVIntervalTree<String> partitions} is not serializable.
    // Please do not change signature (including Object return type) as that may break Serialization.
    private Object writeReplace() {
        return new Serialized(this);
    }

    /**
     * Serialized form of {@link ShardPartitioner} instances.
     * <p>
     *     Solves the issue with non-Java serializable {@link SVIntervalTree}.
     * </p>
     * <p>
     *     For some unknown reason currently we cannot apply Kryo serialization on classes that
     *     extends {@link Serializable} and since {@link ShardPartitioner} extends {@link Partitioner}
     *     we don't have an option but to use Java's standard object serialization framework.
     * </p>
     */
    private static class Serialized implements Serializable {

        private static final long serialVersionUID = ShardPartitioner.serialVersionUID;

        private final String clazzName;
        private final int numberOfPartitions;
        private final String[] contigNames;
        private final int[] contigIndex;
        private final int[] intervalStart;
        private final int[] intervalEnd;
        private final int[] intervalIndex;
        private final int[] intervalPartition;

        private Serialized(final ShardPartitioner<?> partitioner) {
           this.clazzName = partitioner.clazz.getName();
           this.numberOfPartitions = partitioner.numberOfPartitions;
           this.contigNames = new String[partitioner.contigToIndex.size()];
           this.contigIndex = new int[this.contigNames.length];
           int nextIndex = 0;
           for (final Map.Entry<String, Integer> entry : partitioner.contigToIndex.entrySet()) {
               contigNames[nextIndex] = entry.getKey();
               contigIndex[nextIndex++] = entry.getValue();
           }
           this.intervalStart = new int[partitioner.partitions.size()];
           this.intervalEnd = new int[intervalStart.length];
           this.intervalIndex = new int[intervalStart.length];
           this.intervalPartition = new int[intervalStart.length];
           nextIndex = 0;
           for (final SVIntervalTree.Entry<Integer> entry : partitioner.partitions) {
               final SVInterval interval = entry.getInterval();
               intervalStart[nextIndex] = interval.getStart();
               intervalEnd[nextIndex] = interval.getEnd();
               intervalIndex[nextIndex] = interval.getContig();
               intervalPartition[nextIndex++] = entry.getValue();
           }
        }

        // Please don't change signature including return type.
        private Object readResolve() {
            return new ShardPartitioner<>(this);
        }
    }
}
