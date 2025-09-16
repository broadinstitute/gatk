package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;

import java.util.*;

public class SVFederationCollapser extends CanonicalSVCollapser {
    public SVFederationCollapser(ReferenceSequenceFile reference,
                                 AltAlleleSummaryStrategy altAlleleSummaryStrategy,
                                 BreakpointSummaryStrategy breakpointSummaryStrategy,
                                 FlagFieldLogic flagFieldLogic) {
        super(reference, altAlleleSummaryStrategy, breakpointSummaryStrategy, flagFieldLogic);
    }

    @Override
    protected Set<String> collapseFilters(final List<SVCallRecord> items) {
        final Set<String> mergedFilters = new HashSet<>();
        for (final SVCallRecord item : items) {
            final Set<String> filters = item.getFilters();
            // if any member is PASS, set to PASS
            if (filters.isEmpty()) { return Collections.emptySet(); }
            // otherwise return union of all filters
            mergedFilters.addAll(filters);
        }
        return mergedFilters;
    }
}
