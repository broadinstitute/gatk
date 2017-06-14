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
 * a finite false positive probability can be tolerated. Optimal index size and number of hashes are determined
 * for a given number of elements to be inserted and false prositive probability.
 */
@DefaultSerializer(LongBloomFilter.Serializer.class)
public final class LongBloomFilter {

    private final transient Logger logger = LogManager.getLogger(this.getClass());

    private final long totalBits; //Bloom filter bits
    private final int numHashes; //Number of Bloom filter hash functions
    private final long totalBuckets; //Number of 8-bit buckets
    private final int numBucketArrays; //Number of arrays of buckets
    private final int bucketArraySize; //Size of each bucket array (except the last one)
    private final int finalBucketArraySize; //Size of the last bucket array
    private final byte[][] buckets;

    //Allowable number of bits. These are primes near powers of 2. Supports up to ~100GB.
    private final static long[] legalBitSizes = {10007L, 16411L, 32771L, 65537L,
            131101L, 262147L, 524309L, 1048583L, 2097169L,
            4194319L, 8388617L, 16777259L, 33554467L, 67108879L,
            134217757L, 268435459L, 536870923L, 1073741827L, 2147483659L,
            4294967311L, 8589934609L, 17179869209L, 34359738337L, 68719476767L,
            137438953481L, 274877906951L, 549755813881L, 1099511627791L};

    private final static long HASH_SEED_2 = 0x6cebe6dca7f118a6L;

    public LongBloomFilter(final long numElements, final double fpp) {
        Utils.validateArg(numElements > 0, "Number of elements must be greater than 0");
        Utils.validateArg(fpp > 0 && fpp < 1, "False positive probability must be between 0 and 1");

        final long optimalNumberOfBits = getOptimalNumberOfBits(numElements, fpp);

        long smallestLegalSize = -1;
        for (final long legalSize : legalBitSizes) {
            if (legalSize > optimalNumberOfBits) {
                smallestLegalSize = legalSize;
                break;
            }
        }
        if (smallestLegalSize == -1) {
            throw new GATKException("Could not create Bloom filter with " + optimalNumberOfBits + " bits");
        }
        totalBits = smallestLegalSize;

        final int optimalNumberOfHashes = (int) Math.ceil(-Math.log(fpp) / Math.log(2));
        numHashes = optimalNumberOfHashes > 0 ? optimalNumberOfHashes : 1;

        totalBuckets = ((totalBits / 8) + (totalBits % 8 > 0 ? 1 : 0));
        bucketArraySize = SetSizeUtils.legalSizes[SetSizeUtils.legalSizes.length - 1];
        numBucketArrays = (int) (totalBuckets / bucketArraySize) + 1;
        finalBucketArraySize = (int) (totalBuckets % bucketArraySize);

        buckets = new byte[numBucketArrays][];
        for (int i = 0; i < numBucketArrays - 1; i++) {
            buckets[i] = new byte[bucketArraySize];
        }
        buckets[numBucketArrays - 1] = new byte[finalBucketArraySize];
    }

    protected LongBloomFilter(final Kryo kryo, final Input input) {

        totalBits = input.readLong();
        totalBuckets = input.readLong();
        numBucketArrays = input.readInt();
        bucketArraySize = input.readInt();
        finalBucketArraySize = input.readInt();
        numHashes = input.readInt();
        buckets = new byte[numBucketArrays][];
        for (int i = 0; i < numBucketArrays - 1; i++) {
            buckets[i] = input.readBytes(bucketArraySize);
        }
        buckets[numBucketArrays - 1] = input.readBytes(finalBucketArraySize);

        if (logger.isDebugEnabled()) {
            final long X = countBits();
            final long estimatedSize = estimateSize(X);
            logger.debug("Deserialized: totalBits : " + totalBits + ", totalBuckets: " + totalBuckets + ", numHashes: " + numHashes + ", bits set: " + X + ", approximate size: " + estimatedSize);
        }
    }

    protected void serialize(final Kryo kryo, final Output output) {
        output.writeLong(totalBits);
        output.writeLong(totalBuckets);
        output.writeInt(numBucketArrays);
        output.writeInt(bucketArraySize);
        output.writeInt(finalBucketArraySize);
        output.writeInt(numHashes);
        for (int i = 0; i < numBucketArrays; i++) {
            output.writeBytes(buckets[i]);
        }
    }

    public static long getOptimalNumberOfBits(final long numElements, final double fpp) {
        return (long) Math.ceil(-numElements * Math.log(fpp) / (Math.log(2) * Math.log(2)));
    }

    public double getTheoreticalFPP(final long numElements) {
        return Math.pow(1.0 - Math.pow(1.0 - (1.0/totalBits), numHashes * numElements), numHashes);
    }

    private long countBits() {
        final int[] bitCountMap = new int[256];
        for (int b = 0; b < 256; b++) {
            bitCountMap[b] = Integer.bitCount((byte)b);
        }
        long sum = 0;
        for (int i = 0; i < numBucketArrays; i++) {
            for (final byte b : buckets[i]) {
                sum += bitCountMap[0 & b];
            }
        }
        return sum;
    }

    /**
     * Returns an estimate for the number of longs that have been added to the set
     */
    private long estimateSize(final long numBitsSet) {
        return (long) (-(((double) totalBits) / numHashes) * Math.log(1 - (((double) numBitsSet) / totalBits)));
    }

    public boolean add(final long entryValue) {
        final long hash1 = SVUtils.fnvLong64(entryValue);
        final long hash2 = SVUtils.fnvLong64(HASH_SEED_2, entryValue);
        for (int i = 0; i < numHashes; i++) {
            final long bitIndex = applyHashFunction(i, hash1, hash2);
            final int bucketArray = bitIndexToBucketArray(bitIndex);
            final int bucketIndex = bitIndexToBucketIndex(bitIndex);
            buckets[bucketArray][bucketIndex] |= bucketMask(bitIndex);
        }
        return true;
    }

    public boolean contains(final long key) {
        final long hash1 = SVUtils.fnvLong64(key);
        final long hash2 = SVUtils.fnvLong64(HASH_SEED_2, key);
        for (int i = 0; i < numHashes; i++) {
            final long bitIndex = applyHashFunction(i, hash1, hash2);
            final int bucketArray = bitIndexToBucketArray(bitIndex);
            final int bucketIndex = bitIndexToBucketIndex(bitIndex);
            if ((bucketMask(bitIndex) & buckets[bucketArray][bucketIndex]) == 0) return false;
        }
        return true;
    }

    public void addAll(final long[] entryValues) {
        for (final long val : entryValues) {
            add(val);
        }
    }

    public boolean containsAll(final long[] vals) {
        for (final long val : vals) {
            if (!contains(val))
                return false;
        }
        return true;
    }

    /**
     * Computes ith hash. First, two hashes (h_, h_2) are computed using FNV1.
     * Hash i is then computed by double hashing: h_i = h_1 + i * h_2. See reference:
     * <p>
     * Kirsch and Mitzenmacher. 2008. Less hashing, same performance: Building a better Bloom filter. Random
     * Structures & Algorithms. 33:2, 187-218.
     */
    private long applyHashFunction(final int i, final long fnvHash1, final long fnvHash2) {
        final long result = (fnvHash1+ i * fnvHash2) % totalBits;
        return result < 0 ? result + totalBits : result;
    }

    /**
     * Determines partition corresponding to the given the bit index
     */
    private int bitIndexToBucketArray(final long bitIndex) {
        return (int) ((bitIndex >>> 3) / bucketArraySize);
    }

    /**
     * Determines index within the partition corresponding to the given the bit index
     */
    private int bitIndexToBucketIndex(final long bitIndex) {
        return (int) ((bitIndex >>> 3) % bucketArraySize);
    }

    /**
     * Returns bucket bit mask with 1 in the position of given bit index
     */
    private byte bucketMask(final long bitIndex) {
        return (byte) (1 << (bitIndex & 7));
    }

    public void clear() {
        for (int i = 0; i < numBucketArrays; i++) {
            Arrays.fill(buckets[i], (byte) 0);
        }
    }

    public boolean isEmpty() {
        for (final byte[] array : buckets) {
            for (final byte b : array) {
                if (b != 0) return false;
            }
        }
        return true;
    }

    @SuppressWarnings("SimplifiableIfStatement")
    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (!(o instanceof LongBloomFilter)) return false;

        final LongBloomFilter that = (LongBloomFilter) o;

        if (totalBits != that.totalBits) return false;
        if (numHashes != that.numHashes) return false;
        if (numBucketArrays != that.numBucketArrays) return false;
        for (int i = 0; i < numBucketArrays; i++) {
            if (!Arrays.equals(buckets[i], that.buckets[i])) return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int result = (int) (totalBits ^ (totalBits >>> 32));
        result = 31 * result + numHashes;
        for (final byte[] array : buckets) {
            result = 31 * result + Arrays.hashCode(array);
        }
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
