package org.broadinstitute.hellbender.tools.variantdb;

public class BinaryUtils {
    // Function to extract k bits from p position (0-based) 
    // and returns the extracted value as integer 
    static long extractBits(long number, int p, int k) { 
        // make a bit-mask of the desired number of bits
        long mask = ((1L << k) - 1L);
       
        // shift desired data to be the lowest ordered bits, and apply mask
        return (mask) & (number >>> p); 
    } 

    static long encodeTo8Bits(float e, float minValue, float maxValue) {
        float range = maxValue - minValue;
        float n = (e - minValue) / range;
        return Math.round(n * 256.0f);
    } 


    static float decodeFrom8Bits(int i, float minValue, float maxValue) {
        float range = maxValue - minValue;
        float n = (1.0f / 256.0f) * ((float) i);
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
