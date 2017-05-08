package org.broadinstitute.hellbender.tools.spark.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.spark.sv.SVUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Random;

/**
 * Set of long's that is too large to fit in a single LongBloomFilter. Works by binning each entry into a
 * corresponding LongBloomFilter partition using the entry hash value. Partitioning is handled internally given
 * a maximum partition size, number of elements to be added, and false positive probability.
 */
@DefaultSerializer(LargeLongBloomFilter.Serializer.class)
public class LargeLongBloomFilter implements LargeQueryableLongSet {

    private final transient Logger logger = LogManager.getLogger(this.getClass());
    private List<LongBloomFilter> filters;
    private int numFilters;

    public LargeLongBloomFilter(final long maxPartitionBytes, final long numElements, final double fpp) {
        initialize(maxPartitionBytes, numElements, fpp);
    }

    @SuppressWarnings("unchecked")
    protected LargeLongBloomFilter(final Kryo kryo, final Input stream) {
        final boolean oldReferences = kryo.getReferences();
        kryo.setReferences(false);

        numFilters = stream.readInt();
        filters = new ArrayList<>(numFilters);
        for (int i = 0; i < numFilters; i++) {
            filters.add((LongBloomFilter) kryo.readClassAndObject(stream));
        }

        kryo.setReferences(oldReferences);
    }

    /**
     * Our hash must be different than that used in LongBloomFilter
     * Here we use an offset ("start") value that is the FNV-1a hash of "pathseq"
     */
    private static int longHash(final long entryVal) {
        return SVUtils.fnvLong(0x1035e606, entryVal);
    }

    /**
     * Estimates the false positive probability using random values
     */
    public static double estimateFalsePositiveProbability(final LargeLongBloomFilter bf, final long numTrials, final long seed) {
        final Random rd = new Random(seed);
        int numPositive = 0;
        for (int i = 0; i < numTrials; i++) {
            if (bf.contains(rd.nextLong() >>> 1)) {
                numPositive++;
            }
        }
        return numPositive / (double) numTrials;
    }

    private void initialize(final long maxPartitionBytes, final long numElements, final double fpp) {
        Utils.validateArg(maxPartitionBytes > 0, "Max partition size must be at greater than 0");
        Utils.validateArg(numElements > 0, "Number of elements must be greater than 0");
        Utils.validateArg(fpp > 0 && fpp < 1, "False positive probability must be between 0 and 1 (exclusive)");

        logger.debug("Initializing Bloom filter with max partition size " + maxPartitionBytes + " bytes, " + numElements + " elements, and false positive probability " + fpp);

        final long totalBytes = LongBloomFilter.getOptimalNumberOfBytes(numElements, fpp);
        final long legalMaxPartitionBytes = LongBloomFilter.getLegalSizeBelow(maxPartitionBytes);
        final long partitions = (long) Math.ceil(totalBytes / (double) Math.min(totalBytes, legalMaxPartitionBytes));
        final long elementsPerPartition = (long) Math.ceil(numElements / (double) partitions);

        logger.debug("Max legal bytes per partition: " + legalMaxPartitionBytes);
        logger.debug("Optimal number of total bytes: " + totalBytes);
        logger.debug("Optimal number of partitions: " + partitions);
        logger.debug("Elements per partition: " + elementsPerPartition);

        filters = new ArrayList<>((int) partitions);
        for (int i = 0; i < partitions; i++) {
            filters.add(new LongBloomFilter(elementsPerPartition, fpp));
        }
        numFilters = filters.size();
        logger.debug("Created " + numFilters + " partitions");
    }

    protected void serialize(final Kryo kryo, final Output stream) {
        final boolean oldReferences = kryo.getReferences();
        kryo.setReferences(false);

        stream.writeInt(numFilters);
        for (final LongBloomFilter f : filters) {
            kryo.writeClassAndObject(stream, f);
        }

        kryo.setReferences(oldReferences);
    }

    public boolean add(final long entryValue) {
        final int setIndex = setIndexOf(entryValue);
        final LongBloomFilter setRef = filters.get(setIndex);
        return setRef.add(entryValue);
    }

    public void addAll(final long[] entryValues) {
        for (final long val : entryValues) {
            add(val);
        }
    }

    public boolean contains(final long key) {
        return filters.get(setIndexOf(key)).contains(key);
    }

    public boolean containsAll(final long[] vals) {
        for (final long val : vals) {
            if (!contains(val))
                return false;
        }
        return true;
    }

    public void clear() {
        for (final LongBloomFilter f : filters) {
            f.clear();
        }
    }

    public boolean isEmpty() {
        for (final LongBloomFilter s : filters) {
            if (!s.isEmpty()) return false;
        }
        return true;
    }

    public Collection<? extends QueryableLongSet> getSets() {
        return filters;
    }

    /**
     * Computes filter bin for an entry using its hash value
     */
    private int setIndexOf(final long val) {
        return Integer.remainderUnsigned(longHash(val), numFilters);
    }

    @Override
    public int hashCode() {
        int sum = 0;
        for (final LongBloomFilter s : filters) {
            sum += s.hashCode();
        }
        return sum;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final LargeLongBloomFilter that = (LargeLongBloomFilter) o;

        return numFilters == that.numFilters && filters.equals(that.filters);

    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<LargeLongBloomFilter> {
        @Override
        public void write(final Kryo kryo, final Output output, final LargeLongBloomFilter bloomFilter) {
            bloomFilter.serialize(kryo, output);
        }

        @Override
        public LargeLongBloomFilter read(final Kryo kryo, final Input input, final Class<LargeLongBloomFilter> klass) {
            return new LargeLongBloomFilter(kryo, input);
        }
    }

}
