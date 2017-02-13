package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Tests whether a read is "well-formed" -- that is, is free of major internal inconsistencies and issues that could lead
 * to errors downstream. If a read passes this filter, the rest of the hellbender engine should be able to process it without
 * blowing up.
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY)
public final class WellformedReadFilter extends ReadFilter {
    private static final long serialVersionUID = 1l;

    private ReadFilter wellFormedFilter = null;

    // Command line parser requires a no-arg constructor
    public WellformedReadFilter() {
    }

    @Override
    public void setHeader(SAMFileHeader header) {
        super.setHeader(header);
        createFilter();
    }

    public WellformedReadFilter(final SAMFileHeader header) {
        setHeader(header);
    }

    private void createFilter() {
        final AlignmentAgreesWithHeaderReadFilter alignmentAgreesWithHeader = new AlignmentAgreesWithHeaderReadFilter(samHeader);

        wellFormedFilter = ReadFilterLibrary.VALID_ALIGNMENT_START
                .and(ReadFilterLibrary.VALID_ALIGNMENT_END)
                .and(alignmentAgreesWithHeader)
                .and(ReadFilterLibrary.HAS_READ_GROUP)
                .and(ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS)
                .and(ReadFilterLibrary.READLENGTH_EQUALS_CIGARLENGTH)
                .and(ReadFilterLibrary.SEQ_IS_STORED)
                .and(ReadFilterLibrary.CIGAR_CONTAINS_NO_N_OPERATOR);
    }

    @Override
    public boolean test(final GATKRead read ) {
        return wellFormedFilter.test(read);
    }
}
