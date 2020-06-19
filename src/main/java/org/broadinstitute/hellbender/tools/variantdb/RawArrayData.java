package org.broadinstitute.hellbender.tools.variantdb;

import static org.broadinstitute.hellbender.tools.variantdb.BinaryUtils.*;

public class RawArrayData {
    public static enum ArrayGenotype {
        // Order is critical here, the ordinal is the int encoding
        AA,AB, BB, NO_CALL
    }

    // TODO: turn these all into getters/setters with precision checks (e.g. baf)
    int probeId;
    ArrayGenotype genotype;
    Float normx;
    Float normy;
    Float baf;
    Float lrr;

    static ArrayGenotype decodeGenotype(int i) {
        return ArrayGenotype.values()[i];
    }

    static int encodeGenotype(ArrayGenotype g) {
        return g.ordinal();
    }

    public static final int LRR_OFFSET = 0;
    public static final float LRR_MIN = -28;
    public static final float LRR_MAX = 7;

    public static final int BAF_OFFSET = 8;
    public static final float BAF_MIN = 0;
    public static final float BAF_MAX = 1;    

    public static final int NORMX_OFFSET = 16;
    public static final float NORMX_MIN = 0;
    public static final float NORMX_MAX = 8;    

    public static final int NORMY_OFFSET = 24;
    public static final float NORMY_MIN = 0;
    public static final float NORMY_MAX = 8;    

    public static final int GT_OFFSET = 32;
    public static final int PROBE_ID_OFFSET = 42;

    // GTC Data Ranges: https://github.com/Illumina/BeadArrayFiles/blob/develop/docs/GTC_File_Format_v5.pdf
    public static RawArrayData decode(long bits) {

        RawArrayData data = new RawArrayData();
        data.lrr = decodeFrom8Bits((int) extractBits(bits, LRR_OFFSET, 8), LRR_MIN, LRR_MAX);
        data.baf = decodeFrom8Bits((int) extractBits(bits, BAF_OFFSET, 8), BAF_MIN, BAF_MAX);
        data.normx = decodeFrom8Bits((int) extractBits(bits, NORMX_OFFSET, 8), NORMX_MIN, NORMX_MAX);
        data.normy = decodeFrom8Bits((int) extractBits(bits, NORMY_OFFSET, 8), NORMY_MIN, NORMY_MAX);
        data.genotype = decodeGenotype((int) extractBits(bits, GT_OFFSET, 2));
        data.probeId = (int) extractBits(bits, PROBE_ID_OFFSET, 22);

        return data;
    }

    public long encode() {
        long lrrBits = encodeTo8Bits(this.lrr, LRR_MIN, LRR_MAX);
        long bafBits = encodeTo8Bits(this.baf, BAF_MIN, BAF_MAX);
        long normxBits = encodeTo8Bits(this.normx, NORMX_MIN, NORMX_MAX);
        long normyBits = encodeTo8Bits(this.normy, NORMX_MIN, NORMX_MAX);
        long gtBits = (long) encodeGenotype(this.genotype);

        return (
            (lrrBits << LRR_OFFSET) | 
            (bafBits << BAF_OFFSET) |
            (normxBits << NORMX_OFFSET) |
            (normyBits << NORMY_OFFSET) |
            (gtBits << GT_OFFSET) |
            ((long) this.probeId << PROBE_ID_OFFSET )
        );
    }
}
