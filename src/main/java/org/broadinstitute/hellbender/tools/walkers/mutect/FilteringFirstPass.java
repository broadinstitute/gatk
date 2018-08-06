package org.broadinstitute.hellbender.tools.walkers.mutect;


import htsjdk.variant.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.List;

public class FilteringFirstPass {
    final List<FilterResult> filterResults;

    public FilteringFirstPass() {
        filterResults = new ArrayList<>();
    }

    public void add(final FilterResult filterResult, final VariantContext vc) {
        filterResults.add(filterResult);
    }

    public List<FilterResult> getFilterResults() {
        return filterResults;
    }

}
