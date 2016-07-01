package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * CustomReadFilter-derived class used for testing custom read filters.
 */
public class TestCustomReadFilter extends CustomReadFilter {
    private static final long serialVersionUID = 1L;

    private SAMFileHeader samHeader;

    public void setSAMFileHeader(SAMFileHeader header) { this.samHeader = header; }

    public boolean test(GATKRead read) { return true;}
}
