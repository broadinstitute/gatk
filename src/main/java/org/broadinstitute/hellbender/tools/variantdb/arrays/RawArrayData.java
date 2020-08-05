package org.broadinstitute.hellbender.tools.variantdb.arrays;

import java.math.*;
import static org.broadinstitute.hellbender.tools.variantdb.BinaryUtils.*;

public class RawArrayData {
    // TODO: turn these all into getters/setters with precision checks (e.g. baf)
    public Float normx;
    public Float normy;
    public Float baf;
    public Float lrr;

    public RawArrayData(Float normx, Float normy, Float baf, Float lrr) {
        this.normx = normx;
        this.normy = normy;
        this.baf = baf;
        this.lrr = lrr;
    }

    public static final int NORMX_OFFSET = 0;
    public static final int NORMY_OFFSET = 16;
    public static final int LRR_OFFSET = 32;
    public static final int BAF_OFFSET = 48;

    private static final int MIN_16_BIT_VALUE = 0;
    private static final int MAX_16_BIT_VALUE = (int) Math.pow(2, 16) - 2; // reserve for null
    private static final int NULL_ENCODING =  MAX_16_BIT_VALUE + 1;

    private static final int MIN_10_BIT_VALUE = 0;
    private static final int MAX_10_BIT_VALUE = (int) Math.pow(2, 10) - 2; // reserve for null
    private static final int NULL_10_BIT_ENCODING =  MAX_10_BIT_VALUE + 1;

    // store a float with 3-decimal digits in 16 bits by
    // multiplying by 1000 and capping values, reserving 
    // xFFFF FFFF to represent null
    public static int encode(Float f) {
        return encode(f,0);
    }
    public static int encode(Float f, float offset) {

        // TODO: fix to be 10-bit null encoding also...
        if (f == null) return NULL_ENCODING;

        return
            Math.min(
                Math.max(
                    Math.round((f+offset) * 1000.0f),
                    MIN_16_BIT_VALUE
                ),
                MAX_16_BIT_VALUE
            );
    }

    public static int encode10bits(Float f, float offset) {

        if (f == null) return NULL_10_BIT_ENCODING;

        return
            Math.min(
                Math.max(
                    Math.round((f+offset) * 1000.0f),
                    MIN_10_BIT_VALUE
                ),
                MAX_10_BIT_VALUE
            );
    }

    public static Float decode(long i) {
        return decode(i, 0);
    }

    public static Float decode(long i, float offset) {
        if (i == NULL_ENCODING) return null;
        return (
            (float) i) / 1000.0f - offset;
    }

    public static Float decode10bits(long i) {
        if (i == NULL_10_BIT_ENCODING) return null;
        return (
            (float) i) / 1000.0f;
    }

    public RawArrayData(long bits) {
        try {
            this.normx = decode(extractBits(bits, NORMX_OFFSET, 16));
            this.normy = decode(extractBits(bits, NORMY_OFFSET, 16));
            this.lrr = decode(extractBits(bits, LRR_OFFSET, 16), 32.0f);
            this.baf = decode10bits(extractBits(bits, BAF_OFFSET, 10));
        } catch (NullPointerException npe) {
            npe.printStackTrace();
            throw npe;
        }
    }

    public long encode() {
        long normxBits = encode(this.normx);
        long normyBits = encode(this.normy);
        long lrrBits = encode(this.lrr, 32.0f);
        long bafBits = encode10bits(this.baf, 0.0f);

        return (
            (normxBits << NORMX_OFFSET) | 
            (normyBits << NORMY_OFFSET) |
            (lrrBits   << LRR_OFFSET) |
            (bafBits   << BAF_OFFSET) 
        );
    }
}
