package org.broadinstitute.hellbender.tools.variantdb.arrays;

import static org.broadinstitute.hellbender.tools.variantdb.BinaryUtils.*;

import org.broadinstitute.hellbender.exceptions.GATKException;

public class BasicArrayData {
    // TODO remove this before prod
    // replace with
    public static enum ArrayGenotype {
        // Order is critical here, the ordinal is the int encoding
        AA,AB, BB, NO_CALL
    }

    public int sampleId;
    public int probeId;
    public ArrayGenotype genotype;

    public static final int GT_LENGTH = 2;
    public static final int PROBE_ID_LENGTH = 30;
    public static final int MAX_PROBE_ID_VALUE = (int) Math.pow(2, PROBE_ID_LENGTH) - 1;

    public static final int SAMPLE_ID_LENGTH = 30;
    public static final int MAX_SAMPLE_ID_VALUE = (int) Math.pow(2, SAMPLE_ID_LENGTH) - 1;

    public static final int GT_OFFSET = 0;
    public static final int PROBE_ID_OFFSET = GT_OFFSET + GT_LENGTH;
    public static final int SAMPLE_ID_OFFSET = PROBE_ID_OFFSET + PROBE_ID_LENGTH;

    public BasicArrayData(int sampleId, int probeId, ArrayGenotype genotype) {
        // check that the sizes fit
        if (sampleId < 0 || sampleId > MAX_SAMPLE_ID_VALUE) {
            throw new GATKException("Attempted sample id of " + sampleId + " which is great than the maximum of " + MAX_SAMPLE_ID_VALUE);
        }

        if (probeId < 0 || probeId > MAX_PROBE_ID_VALUE) {
            throw new GATKException("Attempted sample id of " + probeId + " which is great than the maximum of " + MAX_PROBE_ID_VALUE);
        }

        this.sampleId = sampleId;
        this.probeId = probeId;
        this.genotype = genotype;
    }

    public BasicArrayData(long bits) {
        this.genotype = decodeGenotype((int) extractBits(bits, GT_OFFSET, GT_LENGTH));
        this.probeId = (int) extractBits(bits, PROBE_ID_OFFSET, PROBE_ID_LENGTH);
        this.sampleId = (int) extractBits(bits, SAMPLE_ID_OFFSET, SAMPLE_ID_LENGTH);
    }

    private ArrayGenotype decodeGenotype(int i) {
        return ArrayGenotype.values()[i];
    }

    private int encodeGenotype(ArrayGenotype g) {
        return g.ordinal();
    }

    public long encode() {
        return (
            ((long) encodeGenotype(this.genotype) << GT_OFFSET) |
            ((long) this.probeId << PROBE_ID_OFFSET ) |
            ((long) this.sampleId << SAMPLE_ID_OFFSET )
        );
    }
}
