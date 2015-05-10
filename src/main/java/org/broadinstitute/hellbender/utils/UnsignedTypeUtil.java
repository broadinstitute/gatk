package org.broadinstitute.hellbender.utils;

/**
 * A utility class for dealing with unsigned types.  This class is primarily used for promoting a value in an unsigned type to
 * the same value in the net larger type of the same form (e.g. Integer to Long)
 */
public final class UnsignedTypeUtil {

    /** Convert an unsigned byte to a signed int */
    public static int uByteToInt(final byte unsignedByte) {
        return unsignedByte & 0xFF;
    }

    /** Convert an unsigned byte to a signed short */
    public static int uByteToShort(final byte unsignedByte) {
        return (short) unsignedByte & 0xFF;
    }

    /** Convert an unsigned short to an Int */
    public static int uShortToInt(final short unsignedShort) {
        return unsignedShort & 0xFFFF;
    }

    /** Convert an unsigned int to a long */
    public static long uIntToLong(final int unsignedInt) {
        return unsignedInt & 0xFFFFFFFFL;
    }
}
