package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.util.Locatable;

/**
 * Indicates that this object has a genomic location and provides a systematic interface to get it.
 */
public interface HasGenomeLocation {
    public Locatable getLocation();
}
