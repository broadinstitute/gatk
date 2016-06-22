package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;

/**
 * Interface to be implemented by custom read filters that can be dynamically
 * discovered at runtime and instantiated by using the "CUSTOM_READ_FILTER"
 * CommandLineReadFilter.
 */
public abstract class CustomReadFilter implements ReadFilter {

    private static final long serialVersionUID = 1L;

    /**
     * Set the SAMFileHeader being used for this read filter.
     * @param header
     */
    abstract public void setSAMFileHeader(SAMFileHeader header);
}
