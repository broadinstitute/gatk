package org.broadinstitute.hellbender.utils;

/**
 * Indicates that this object has a genomic location and provides a systematic interface to get it.
 */
public interface HasGenomeLocation {
    public GenomeLoc getLocation();
}
