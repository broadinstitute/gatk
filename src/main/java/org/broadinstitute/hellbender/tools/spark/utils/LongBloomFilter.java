package org.broadinstitute.hellbender.tools.spark.utils;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.SVUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;

/**
 * Bloom filter for primitive longs. Useful for quickly querying whether a long is a part of a set when a
 * a finite false positive probability can be tolerated. Optimal array size and number of hashes are determined
 * for a given number of elements to be inserted and false prositive probability.
 */
@DefaultSerializer(LongBloomFilter.Serializer.class)
public class LongBloomFilter implements QueryableLongSet {

    // largest prime numbers less than each half power of 2 from 2^8 to 2^31
    // Note the last size is the greatest prime that is an allowable Java array size (<=2^31-3)
    private static final int[] legalSizes = {
            251, 359, 509, 719, 1021, 1447, 2039, 2887, 4093, 5791, 8191, 11579, 16381, 23167, 32749, 46337, 65521,
            92681, 131071, 185363, 262139, 370723, 524287, 741431, 1048573, 1482907, 2097143, 2965819, 4194301, 5931641,
            8388593, 11863279, 16777213, 23726561, 33554393, 47453111, 67108859, 94906249, 134217689, 189812507,
            268435399, 379625047, 536870909, 759250111, 1073741789, 1518500213, 2147483629
    };
    private final transient Logger logger = LogManager.getLogger(this.getClass());
    private final long numBits;
    private final int numBuckets, numHashes;
    private final byte[] buckets;

    public LongBloomFilter(final long num_elements, final double fpp) {
        Utils.validateArg(num_elements > 0, "Number of elements must be greater than 0");
        Utils.validateArg(fpp > 0 && fpp <= 1, "False positive probability must be between 0 and 1");

        final long optimal_num_bytes = getOptimalNumberOfBytes(num_elements, fpp);
        this.numBuckets = getLegalSizeAbove(optimal_num_bytes);
        this.numBits = (this.numBuckets * 8L);
        final int optimal_num_hashes = (int) Math.round(-Math.log(fpp) / Math.log(2));
        this.numHashes = optimal_num_hashes > 0 ? optimal_num_hashes : 1;
        this.buckets = new byte[numBuckets];
    }

    public LongBloomFilter(final long minNumElements, final int numBytes) {
        final int maxLegalSize = legalSizes[legalSizes.length-1];
        Utils.validateArg(numBytes > 0 && numBytes <= maxLegalSize, "Number of bytes must be between 0 and " + maxLegalSize);

        this.numBuckets = getLegalSizeAbove(minNumElements);
        this.numBits = (this.numBuckets * 8L);
        final int optimal_num_hashes = (int) Math.ceil(Math.log(2) * minNumElements / this.numBuckets);
        this.numHashes = optimal_num_hashes > 0 ? optimal_num_hashes : 1;
        this.buckets = new byte[numBuckets];
    }

    protected LongBloomFilter(final Kryo kryo, final Input input) {
        final boolean oldReferences = kryo.getReferences();
        kryo.setReferences(false);

        numBits = input.readLong();
        numBuckets = input.readInt();
        numHashes = input.readInt();
        buckets = input.readBytes(numBuckets);

        if (logger.isDebugEnabled()) {
            final long X = countBits();
            final long estimated_size = estimateSize(X);
            logger.debug("Deserialized: numBits : " + numBits + ", numBuckets: " + numBuckets + ", numHashes: " + numHashes + ", bits set: " + X + ", approximate size: " + estimated_size);
        }
        kryo.setReferences(oldReferences);
    }

    private static int getLegalSizeAbove(final long minNumBytes) {
        for (int i = 0; i < legalSizes.length; i++) {
            if (minNumBytes <= legalSizes[i]) {
                return legalSizes[i];
            }
        }
        throw new IllegalArgumentException("No legal sizes large enough for size " + minNumBytes);
    }

    public static int getLegalSizeBelow(final long maxNumBytes) {
        if (maxNumBytes <= legalSizes[0]) {
            throw new GATKException("No legal sizes smaller than size " + maxNumBytes);
        }
        for (int i = 1; i < legalSizes.length; i++) {
            if (maxNumBytes < legalSizes[i]) {
                return legalSizes[i - 1];
            }
        }
        return legalSizes[legalSizes.length - 1];
    }

    public static long getOptimalNumberOfBytes(final long num_elements, final double fpp) {
        return (long) Math.ceil(-num_elements * Math.log(fpp) / (8 * Math.log(2) * Math.log(2)));
    }

    private long countBits() {
        long X = 0;
        for (final byte b : buckets) {
            byte bb = b;
            for (int i = 0; i < Byte.SIZE; i++) {
                X += bb & 1;
                bb >>= 1;
            }
        }
        return X;
    }

    private long estimateSize(final long numBitsSet) {
        return (long) (-(((double) numBits) / numHashes) * Math.log(1 - (((double) numBitsSet) / numBits)));
    }

    protected void serialize(final Kryo kryo, final Output output) {
        final boolean oldReferences = kryo.getReferences();
        kryo.setReferences(false);

        output.writeLong(numBits);
        output.writeInt(numBuckets);
        output.writeInt(numHashes);
        output.writeBytes(buckets);

        kryo.setReferences(oldReferences);
    }

    public boolean add(final long entryValue) {
        final long[] bitIndex = entryIndices(entryValue);
        for (int i = 0; i < bitIndex.length; i++) {
            final int bucket = indexBucket(bitIndex[i]);
            buckets[bucket] |= bucketMask(bitIndex[i]);
        }
        return true;
    }

    public void addAll(final long[] entryValues) {
        for (final long val : entryValues) {
            add(val);
        }
    }

    /**
     * Using different semantics than a Set because false positives are possible
     */
    public boolean contains(final long key) {
        final long[] bitIndex = entryIndices(key);
        for (int i = 0; i < bitIndex.length; i++) {
            final int bucket = indexBucket(bitIndex[i]);
            if ((bucketMask(bitIndex[i]) & buckets[bucket]) == 0) return false;
        }
        return true;
    }

    public boolean containsAll(final long[] vals) {
        for (final long val : vals) {
            if (!contains(val))
                return false;
        }
        return true;
    }

    /**
     * Determines partition corresponding to the given the bit index
     */
    private int indexBucket(final long bitIndex) {
        return (int) (bitIndex >>> 3);
    }

    /**
     * Returns bucket bit mask with 1 in the position of given bit index
     */
    private byte bucketMask(final long bitIndex) {
        return (byte) (1 << (bitIndex & 7));
    }

    public void clear() {
        for (int i = 0; i < numBuckets; i++) {
            buckets[i] = 0;
        }
    }

    public boolean isEmpty() {
        for (final byte b : buckets) {
            if (b != 0) return false;
        }
        return true;
    }

    /**
     * Computes hash values for a given entry. First, two hashes (h_, h_2) are computed using FNV1.
     * Hash i is then computed by double hashing: h_i = h_1 + i * h_2. See reference:
     * <p>
     * Kirsch and Mitzenmacher. 2008. Less hashing, same performance: Building a better Bloom filter. Random
     * Structures & Algorithms. 33:2, 187-218.
     */
    private long[] entryIndices(final long val) {
        final long[] indices = new long[numHashes];
        final long hash = SVUtils.fnvLong64(val);
        final long hash1 = hash >> 32;
        final long hash2 = hash & 0xffffffffL;
        for (int i = 0; i < numHashes; i++) {
            indices[i] = (hash1 + i * hash2) % numBits;
            if (indices[i] < 0) indices[i] += numBits;
        }
        return indices;
    }

    @SuppressWarnings("SimplifiableIfStatement")
    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final LongBloomFilter that = (LongBloomFilter) o;

        if (numBits != that.numBits) return false;
        if (numBuckets != that.numBuckets) return false;
        if (numHashes != that.numHashes) return false;
        return Arrays.equals(buckets, that.buckets);

    }

    @Override
    public int hashCode() {
        int result = (int) (numBits ^ (numBits >>> 32));
        result = 31 * result + numBuckets;
        result = 31 * result + numHashes;
        result = 31 * result + Arrays.hashCode(buckets);
        return result;
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<LongBloomFilter> {
        @Override
        public void write(final Kryo kryo, final Output output, final LongBloomFilter bloomFilter) {
            bloomFilter.serialize(kryo, output);
        }

        @Override
        public LongBloomFilter read(final Kryo kryo, final Input input, final Class<LongBloomFilter> klass) {
            return new LongBloomFilter(kryo, input);
        }
    }
}
