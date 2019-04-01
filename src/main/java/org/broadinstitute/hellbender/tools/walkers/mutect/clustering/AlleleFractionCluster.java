package org.broadinstitute.hellbender.tools.walkers.mutect.clustering;

import java.util.List;

public interface AlleleFractionCluster {
    double logLikelihood(final Datum datum);

    double logLikelihood(final int totalCount, final int altCount);

    void learn(final List<Datum> data);

    @Override
    String toString();
}
