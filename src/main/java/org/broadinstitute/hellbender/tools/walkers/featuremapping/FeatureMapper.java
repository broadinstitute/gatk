package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.function.Consumer;

public interface FeatureMapper {

    void        forEachOnRead(GATKRead read, ReferenceContext referenceContext, Consumer<? super FlowFeatureMapper.MappedFeature> action);
    boolean     noFeatureButFilterAt(GATKRead read, ReferenceContext referenceContext, int start);
}
