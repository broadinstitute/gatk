package org.broadinstitute.hellbender.tools.walkers.mutect.clustering;

import java.util.List;

public interface AlleleFractionCluster {
    double log10Likelihood(final Datum datum);

    double log10Likelihood(final int totalCount, final int altCount);

    void learn(final List<Datum> data);

    @Override
    String toString();
}
