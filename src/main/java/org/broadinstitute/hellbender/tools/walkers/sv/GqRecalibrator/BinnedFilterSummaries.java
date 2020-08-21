package org.broadinstitute.hellbender.tools.walkers.sv.GqRecalibrator;

import java.util.ArrayList;
import java.util.Collection;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

class BinnedFilterSummaries extends ArrayList<FilterSummary> {
    private static final long serialVersionUID = 0;

    BinnedFilterSummaries(final Collection<FilterSummary> filterSummaries) {
        super(filterSummaries);
    }

    BinnedFilterSummaries(final FilterSummary filterSummary, final int binIndex, final int numBins) {
        this(
                IntStream.range(0, numBins)
                        .mapToObj(i -> i == binIndex ? filterSummary : FilterSummary.EMPTY)
                        .collect(Collectors.toList())
        );
    }

    BinnedFilterSummaries add(final BinnedFilterSummaries other) {
        return this.equals(BinnedFilterSummaries.EMPTY) ?
                other :
                other.equals(BinnedFilterSummaries.EMPTY) ?
                        this :
                        new BinnedFilterSummaries(
                                IntStream.range(0, this.size()).mapToObj(
                                        i -> this.get(i).add(other.get(i))).collect(Collectors.toList()
                                )
                        );
    }

    static BinnedFilterSummaries EMPTY = new BinnedFilterSummaries(FilterSummary.EMPTY, 0, 1);
}
