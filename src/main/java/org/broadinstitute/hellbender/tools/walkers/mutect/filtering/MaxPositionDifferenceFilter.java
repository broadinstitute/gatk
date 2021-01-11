package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;

public class MaxPositionDifferenceFilter extends HardFilter {


    public MaxPositionDifferenceFilter() {
    }

    @Override
    public ErrorType errorType() { return ErrorType.ARTIFACT; }

    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringEngine filteringEngine) {
        final List<Integer> readPositionMaxDiff = vc.getAttributeAsIntList(GATKVCFConstants.READ_POSITION_MAX_DIFF_KEY, 0);
        final double[] tumorLods = Mutect2FilteringEngine.getTumorLogOdds(vc);
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);

        return readPositionMaxDiff.get(indexOfMaxTumorLod + 1) == 1;
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.MAX_POSITION_DIFFERENCE_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() { return Collections.singletonList(GATKVCFConstants.READ_POSITION_MAX_DIFF_KEY); }
}
