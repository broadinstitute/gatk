package org.broadinstitute.hellbender.engine.filters.flow;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * A read filter to test if the TP values for each hmer in a flow based read form
 * are wihin the allowed range (being the possible lengths of hmers - maxHmer)
 */public class FlowBasedTPAttributeValidReadFilter extends ReadFilter implements FlowBasedHmerBasedReadFilterHelper.FilterImpl {
    private static final long serialVersionUID = 1l;

    @Argument(fullName = "read-filter-max-hmer",
            doc = "maxHmer to use for testing in the filter", optional = true)
    public int maxHmer = 12;

    public FlowBasedTPAttributeValidReadFilter() {
        super();
    }

    @Override
    public boolean test(final GATKRead read) {
        return FlowBasedHmerBasedReadFilterHelper.test(read, this);
    }

    @Override
    public byte[] getValuesOfInterest(final GATKRead read) {
        return read.getAttributeAsByteArray(FlowBasedRead.FLOW_MATRIX_TAG_NAME);
    }

    @Override
    public boolean testHmer(final byte[] values, final int hmerStartingOffset, final int hmerLength) {

        // check matrix index resulting from tp value does not exceed limits
        // (note that tp value is a 1/0/-1 adjustment of the hmer length
        for ( int i = 0 ; i < hmerLength ; i++ ) {
            final int targetValue = values[hmerStartingOffset + i] + hmerLength;

            if (targetValue < 0 || targetValue > maxHmer)
                return false;
        }

        return true;
    }
}
