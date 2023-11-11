package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.function.Consumer;

public interface FeatureMapper {

    enum FilterStatus {
        None,
        Filtered,
        NoFeatureAndFiltered
    };

    void        forEachOnRead(GATKRead read, ReferenceContext referenceContext, Consumer<? super FlowFeatureMapper.MappedFeature> action);
    FilterStatus noFeatureButFilterAt(GATKRead read, ReferenceContext referenceContext, int start);
}
