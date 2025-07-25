package org.broadinstitute.hellbender.tools.sv.concordance;

import org.broadinstitute.hellbender.tools.sv.SVCallRecord;

import java.util.List;

public interface SVMatcher {
    List<ClosestSVFinder.LinkageConcordanceRecord> flush(final boolean force);

    String getLastItemContig();

    void add(final SVCallRecord item, final Long id, final boolean isTruthVariant);


}
