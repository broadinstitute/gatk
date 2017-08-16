package org.broadinstitute.hellbender.utils.iterators;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * Some utility routines for arrays.
 */
public class ArrayUtils {

    /**
     * Returns an array of ints with the same values as the input array of bytes.
     * @param bytes the input bytes.
     * @param signed whether the input bytes should treated as signed or unsigned
     *               values.
     * @return never {@code null}.
     */
    public static int[] toInts(final byte[] bytes, final boolean signed) {
        Utils.nonNull(bytes);
        final int[] result = new int[bytes.length];
        if (signed) {
            for (int i = 0; i < result.length; i++) {
                result[i] = bytes[i];
            }
        } else {
            for (int i = 0; i < result.length; i++) {
                result[i] = 0x000000FF & bytes[i];
            }
        }
        return result;
    }
}
