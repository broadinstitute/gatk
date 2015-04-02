package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * Checks to ensure that the alignment of each read makes sense based on the contents of the header.
 */
public final class AlignmentAgreesWithHeaderReadFilter implements ReadFilter {
    private static final long serialVersionUID = 1l;

    private final SAMFileHeader header;

    public AlignmentAgreesWithHeaderReadFilter( final SAMFileHeader header ) {
        this.header = header;
    }

    @Override
    public boolean test( GATKRead read ) {
        return ReadUtils.alignmentAgreesWithHeader(header, read);
    }
}
