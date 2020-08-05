package org.broadinstitute.hellbender.tools.variantdb;

public class BinaryUtils {
    // Function to extract k bits from p position (0-based) 
    // and returns the extracted value as integer 
    public static long extractBits(long number, int p, int k) {
        // make a bit-mask of the desired number of bits
        long mask = ((1L << k) - 1L);
       
        // shift desired data to be the lowest ordered bits, and apply mask
        return (mask) & (number >>> p); 
    } 

    // 0xFF (255) is reserved as NULL
    public static long encodeTo8Bits(Float e, float minValue, float maxValue) {
        if (e == null) {
            return 255;
        }

        if (e > maxValue) {
            e = maxValue;
        }

        if (e < minValue) {
            e = minValue;
        }

        float range = maxValue - minValue;
        float n = (e - minValue) / range;
        return Math.round(n * 254.0f);
    } 


    // 0xFF (255) is reserved as NULL
    public static Float decodeFrom8Bits(int i, float minValue, float maxValue) {
        if (i == 255) {
            return null;
        }

        float range = maxValue - minValue;
        float n = (1.0f / 254.0f) * ((float) i);
        return n * range + minValue;
    } 

    /**
     * Converts an long to a 64-bit binary string
     * @param number
     *      The number to convert
     * @param groupSize
     *      The number of bits in a group
     * @return
     *      The 64-bit long bit string
     */
    public static String longToBinaryString(long number, int groupSize) {
        StringBuilder result = new StringBuilder();

        for(long i = 63; i >= 0 ; i--) {
            long mask = 1L << i;
            result.append((number & mask) != 0 ? "1" : "0");

            if (i % groupSize == 0)
                result.append(" ");
        }
        result.replace(result.length() - 1, result.length(), "");

        return result.toString();
    }
}
