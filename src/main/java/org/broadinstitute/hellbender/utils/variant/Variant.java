package org.broadinstitute.hellbender.utils.variant;

import htsjdk.samtools.util.Locatable;

import java.util.UUID;

/**
 * Variant is (currently) a minimal variant interface needed by the Hellbender pipeline.
 * This will be expanded as more methods are needed.
 * NOTE: getStart() and getEnd() are 1-base inclusive (which matches the current GATK tools).
 * This does not match the GA4GH spec.
 */
public interface Variant extends Locatable {
    boolean isSnp();
    boolean isIndel();
    UUID getUUID();
}
