package org.broadinstitute.hellbender.engine.filters.flow;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * A read filter to test if the read's readGroup has a flow order associated with it
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out reads without a flow associated")
public class ReadGroupHasFlowOrderReadFilter extends ReadFilter {
    private static final long serialVersionUID = 1l;
    private final static OneShotLogger readGroupFiltered = new OneShotLogger(ReadGroupHasFlowOrderReadFilter.class);

    public ReadGroupHasFlowOrderReadFilter() {

    }

    public ReadGroupHasFlowOrderReadFilter(final SAMFileHeader header) {
        setHeader(header);
    }

    @Override
    public boolean test(final GATKRead read) {

        final boolean     result;

        if ( read.getReadGroup() == null ) {
            result = false;
        } else if ( samHeader.getReadGroup(read.getReadGroup()) == null ) {
            result = false;
        } else if ( samHeader.getReadGroup(read.getReadGroup()).getFlowOrder() == null ) {
            result = false;
        } else {
            result = true;
        }

        if ( !result ) {
            readGroupFiltered.warn("at least one of  readgroup is missing or missing a flow order.");
        }

        return result;
    }
}
