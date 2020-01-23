package org.broadinstitute.hellbender.tools.walkers.mutect.clustering;

import java.util.List;

public interface AlleleFractionCluster {
    /**
     * The log likelihood of real variation relative to a sequencing error log likelihood of zero
     * obtained by correcting the TLOD of Mutect2's somatic likelihoods model, which has a flat prior on allele fraction,
     * to account for a clustered allele fraction.
     */
    double correctedLogLikelihood(final Datum datum);

    double logLikelihood(final int totalCount, final int altCount);

    void learn(final List<Datum> data, final double[] responsibilities);

    @Override
    String toString();
}
