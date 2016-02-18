package org.broadinstitute.hellbender.utils.variant;

import htsjdk.samtools.util.Locatable;

/**
 * Variant is (currently) a minimal variant interface needed by the Hellbender pipeline.
 * This will be expanded as more methods are needed.
 * NOTE: getStart() and getEnd() are 1-base inclusive (which matches the current GATK tools).
 * This does not match the GA4GH spec.
 */
public interface GATKVariant extends Locatable {
    boolean isSnp();
    boolean isIndel();
}
