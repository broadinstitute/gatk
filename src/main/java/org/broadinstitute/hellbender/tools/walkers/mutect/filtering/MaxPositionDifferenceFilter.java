package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MaxPositionDifferenceFilter extends HardFilter {


    public MaxPositionDifferenceFilter() {
    }

    @Override
    public ErrorType errorType() {
        return ErrorType.ARTIFACT;
    }

    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringEngine filteringEngine) {
        final List<Integer> startPositionMaxDiff = vc.getAttributeAsIntList(GATKVCFConstants.READ_START_POSITION_MAX_DIFF_KEY, 0);
        final List<Integer> startPositionMinDiff = vc.getAttributeAsIntList(GATKVCFConstants.READ_START_POSITION_MIN_DIFF_KEY, 0);

        final double[] tumorLods = Mutect2FilteringEngine.getTumorLogOdds(vc);
        final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);

        final int threshold = 3;

        if (startPositionMaxDiff.get(indexOfMaxTumorLod) <= threshold && startPositionMinDiff.get(indexOfMaxTumorLod) <= threshold) {
            return true;
        } else {
            return false;
        }
    }

    @Override
    public String filterName() {
        return GATKVCFConstants.MAX_POSITION_DIFFERENCE_FILTER_NAME;
    }

    @Override
    protected List<String> requiredAnnotations() {
        return new ArrayList<>(Arrays.asList(GATKVCFConstants.READ_START_POSITION_MAX_DIFF_KEY,
                                             GATKVCFConstants.READ_START_POSITION_MIN_DIFF_KEY,
                                             GATKVCFConstants.INSERT_SIZE_DIFF_KEY));
    }
}