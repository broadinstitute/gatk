package org.broadinstitute.hellbender.engine.filters.flow;

import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * A read filter to test if the TP values for each hmer in a flow based read form
 * a polindrome (as they should)
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out flow reads that have invalid TP values (non palandromic)")
public class FlowBasedTPAttributeSymetricReadFilter extends ReadFilter implements FlowBasedHmerBasedReadFilterHelper.FilterImpl {
    private static final long serialVersionUID = 1l;

    public FlowBasedTPAttributeSymetricReadFilter() {
        super();
    }

    @Override
    public boolean test(final GATKRead read) {

        return FlowBasedHmerBasedReadFilterHelper.test(read, this);
    }

    @Override
    public byte[] getValuesOfInterest(GATKRead read) {
        return read.getAttributeAsByteArray(FlowBasedRead.FLOW_MATRIX_TAG_NAME);
    }

    @Override
    public boolean testHmer(byte[] values, int hmerStartingOffset, int hmerLength) {
        return FlowBasedHmerBasedReadFilterHelper.isPalindrome(values, hmerStartingOffset, hmerLength);
    }

}
