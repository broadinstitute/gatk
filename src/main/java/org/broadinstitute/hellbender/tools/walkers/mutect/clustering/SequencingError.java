package org.broadinstitute.hellbender.tools.walkers.mutect.clustering;

import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.List;

public class SequencingError implements AlleleFractionCluster {

    public SequencingError() { }

    @Override
    public double log10Likelihood(final Datum datum) {
        return 0;
    }

    @Override
    public double log10Likelihood(final int totalCount, final int altCount) {
        throw new GATKException.ShouldNeverReachHereException("This method should never be called on the sequencing error cluster.");
    }

    @Override
    public void learn(final List<Datum> data) { }

    @Override
    public String toString() {
        return "sequencing error";
    }
}
