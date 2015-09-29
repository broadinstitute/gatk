package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Tests whether a read is "well-formed" -- that is, is free of major internal inconsistencies and issues that could lead
 * to errors downstream. If a read passes this filter, the rest of the hellbender engine should be able to process it without
 * blowing up.
 */
public final class WellformedReadFilter implements ReadFilter {
    private static final long serialVersionUID = 1l;

    private final ReadFilter wellFormedFilter;

    public WellformedReadFilter( final SAMFileHeader header ) {
        final AlignmentAgreesWithHeaderReadFilter alignmentAgreesWithHeader = new AlignmentAgreesWithHeaderReadFilter(header);

        wellFormedFilter = ReadFilterLibrary.VALID_ALIGNMENT_START
                             .and(ReadFilterLibrary.VALID_ALIGNMENT_END)
                             .and(alignmentAgreesWithHeader)
                             .and(ReadFilterLibrary.HAS_READ_GROUP)
                             .and(ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS)
                             .and(ReadFilterLibrary.SEQ_IS_STORED)
                             .and(ReadFilterLibrary.CIGAR_IS_SUPPORTED);
    }


    @Override
    public boolean test( final GATKRead read ) {
        return wellFormedFilter.test(read);
    }
}
