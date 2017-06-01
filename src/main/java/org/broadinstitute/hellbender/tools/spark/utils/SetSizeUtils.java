package org.broadinstitute.hellbender.tools.spark.utils;

/**
 * Set size utility
 */
public final class SetSizeUtils {

    /**
     * Largest prime numbers less than each half power of 2 from 2^8 to 2^31
     * Note the last size is the greatest prime that is an allowable Java array size (<=2^31-3)
     */
    public static final int[] legalSizes = {
            251, 359, 509, 719, 1021, 1447, 2039, 2887, 4093, 5791, 8191, 11579, 16381, 23167, 32749, 46337, 65521,
            92681, 131071, 185363, 262139, 370723, 524287, 741431, 1048573, 1482907, 2097143, 2965819, 4194301, 5931641,
            8388593, 11863279, 16777213, 23726561, 33554393, 47453111, 67108859, 94906249, 134217689, 189812507,
            268435399, 379625047, 536870909, 759250111, 1073741789, 1518500213, 2147483629
    };

    /**
     * Computes next largest legal table size for a given number of elements, accounting for load factor
     */
    public static int getLegalSizeAbove(final long minElements, final double loadFactor) {
        final long augmentedSize = (long) (minElements / loadFactor);
        for (final int legalSize : legalSizes) {
            if (legalSize > augmentedSize) return legalSize;
        }
        throw new IllegalArgumentException("No legal sizes large enough for size " + minElements);
    }

    public static int getLegalSizeAbove(final long minElements) {
        return getLegalSizeAbove(minElements, 1.0);
    }

    /**
     * Computes next smallest legal table size for a given number of elements, accounting for load factor
     */
    public static int getLegalSizeBelow(final long maxElements, final double loadFactor) {
        final long augmentedSize = (long) (maxElements / loadFactor);
        if (augmentedSize <= legalSizes[0]) {
            throw new IllegalArgumentException("No legal sizes small enough for size " + maxElements);
        }
        for (int i = 1; i < legalSizes.length; i++) {
            if (augmentedSize <= legalSizes[i]) {
                return legalSizes[i - 1];
            }
        }
        return legalSizes[legalSizes.length - 1];
    }

    public static int getLegalSizeBelow(final long maxElements) {
        return getLegalSizeBelow(maxElements, 1.0);
    }

}
