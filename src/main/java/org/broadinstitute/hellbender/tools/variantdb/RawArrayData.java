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
    float normx;
    float normy;
    float baf;
    float lrr;

    static ArrayGenotype decodeGenotype(int i) {
        return ArrayGenotype.values()[i];
    }

    static int encodeGenotype(ArrayGenotype g) {
        return g.ordinal();
    }

    private static int LRR_OFFSET = 0;
    private static int BAF_OFFSET = 8;
    private static int NORMX_OFFSET = 16;
    private static int NORMY_OFFSET = 24;
    private static int GT_OFFSET = 32;
    private static int PROBE_ID_OFFSET = 42;

    // GTC Data Ranges: https://github.com/Illumina/BeadArrayFiles/blob/develop/docs/GTC_File_Format_v5.pdf
    public static RawArrayData decode(long bits) {

        RawArrayData data = new RawArrayData();
        data.lrr = decodeFrom8Bits((int) extractBits(bits, LRR_OFFSET, 8), -25, 25);
        data.baf = decodeFrom8Bits((int) extractBits(bits, BAF_OFFSET, 8), 0, 1);
        data.normx = decodeFrom8Bits((int) extractBits(bits, NORMX_OFFSET, 8), 0, 65535);
        data.normy = decodeFrom8Bits((int) extractBits(bits, NORMY_OFFSET, 8), 0, 65535);
        data.genotype = decodeGenotype((int) extractBits(bits, GT_OFFSET, 2));
        data.probeId = (int) extractBits(bits, PROBE_ID_OFFSET, 22);

        return data;
    }

    public long encode() {
        return (
            (encodeTo8Bits(this.lrr, -25, 25) << LRR_OFFSET) | 
            (encodeTo8Bits(this.baf, 0, 1) << BAF_OFFSET) |
            (encodeTo8Bits(this.normx, 0, 65535) << NORMX_OFFSET) |
            (encodeTo8Bits(this.normy, 0, 65535) << NORMY_OFFSET) |
            ((long) encodeGenotype(this.genotype) << GT_OFFSET) |
            ((long) this.probeId << PROBE_ID_OFFSET )
        );
    }
}
