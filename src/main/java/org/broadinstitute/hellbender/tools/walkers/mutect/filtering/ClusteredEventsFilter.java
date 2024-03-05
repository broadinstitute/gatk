package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;

public class ClusteredEventsFilter extends HardFilter {
    private final int maxEventsInRegion;

    public ClusteredEventsFilter(final int maxEventsInRegion) {
        this.maxEventsInRegion = maxEventsInRegion;
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringEngine filteringEngine) {
        final List<Integer> eventCounts = vc.getAttributeAsIntList(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY, 0);
        return eventCounts.stream().mapToInt(n -> n).max().getAsInt() > maxEventsInRegion;
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.CLUSTERED_EVENTS_FILTER_NAME;
    }

    @Override
    protected List<String> requiredInfoAnnotations() { return Collections.singletonList(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY); }
}
