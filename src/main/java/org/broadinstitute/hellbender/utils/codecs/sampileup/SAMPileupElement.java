package org.broadinstitute.hellbender.utils.codecs.sampileup;

/**
 * Simple representation of a single base with associated quality from a SAM pileup
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class SAMPileupElement {

    /**
     * The first base
     */
    private final byte base;


    /**
     * Base quality for the first base
     */
    private final byte baseQuality;

    SAMPileupElement(final byte base, final byte baseQuality) {
        this.base = base;
        this.baseQuality = baseQuality;
    }

    /**
     * Get the base
     */
    public byte getBase() {
        return base;
    }

    /**
     * Get the quality
     */
    public byte getBaseQuality() {
        return baseQuality;
    }

}
