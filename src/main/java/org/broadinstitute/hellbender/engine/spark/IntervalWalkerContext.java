package org.broadinstitute.hellbender.engine.spark;

import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;

/**
 * Encapsulates a {@link SimpleInterval} with the reads that overlap it (the {@link ReadsContext} and
 * its {@link ReferenceContext} and {@link FeatureContext}.
 */
public class IntervalWalkerContext {
    private final SimpleInterval interval;
    private final ReadsContext readsContext;
    private final ReferenceContext referenceContext;
    private final FeatureContext featureContext;

    public IntervalWalkerContext(SimpleInterval interval, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        this.interval = interval;
        this.readsContext = readsContext;
        this.referenceContext = referenceContext;
        this.featureContext = featureContext;
    }

    public SimpleInterval getInterval() {
        return interval;
    }

    public ReadsContext getReadsContext() {
        return readsContext;
    }

    public ReferenceContext getReferenceContext() {
        return referenceContext;
    }

    public FeatureContext getFeatureContext() {
        return featureContext;
    }
}
