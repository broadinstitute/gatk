package org.broadinstitute.hellbender.tools.gvs.common;

import java.util.BitSet;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.broadinstitute.hellbender.exceptions.GATKException;

/**
 * Wrapper around a BitSet to track locations where a variant is present in a range of locations.
 * Constructed with a min and max location, so a range and offset can be calculated
 */
public class VariantBitSet {
    private static final Logger logger = LogManager.getLogger(VariantBitSet.class);

    private final BitSet bitSet;
    private final long locationOffset;

    public VariantBitSet(long minLocation, long maxLocation) {
        if (!SchemaUtils.decodeContig(minLocation).equals(SchemaUtils.decodeContig(maxLocation))) {
            throw new GATKException("Can not process cross-contig boundaries");
        }

        // set up a bitset to track locations where there is a variant
        long range = (maxLocation + IngestConstants.MAX_REFERENCE_BLOCK_BASES) - (minLocation - IngestConstants.MAX_REFERENCE_BLOCK_BASES) + 1;
        this.locationOffset = minLocation - Math.max(IngestConstants.MAX_REFERENCE_BLOCK_BASES, IngestConstants.MAX_DELETION_SIZE);

        if (range > Integer.MAX_VALUE) {
            throw new GATKException("Single contig can not be bigger than " + Integer.MAX_VALUE);
        }

        this.bitSet = new BitSet((int) range);
        logger.info("Constructed Variant Filter BitSet using " + range + " bits");
    }

    public void setVariant(long location) {
        bitSet.set((int) (location - locationOffset));
    }

    /**
     * @return true if there was a variant recorded in the range of startLocation (inclusive) to endLocation (exclusive).
     */
    public boolean containsVariant(long startLocation, long endLocation) {
        return !bitSet.get( (int) (startLocation - locationOffset), (int) (endLocation - locationOffset) ).isEmpty();
    }

    // access to underlying bitset and location offset
    public BitSet getBitSet() { return this.bitSet; }
    public long getLocationOffset() { return this.locationOffset; }

}
