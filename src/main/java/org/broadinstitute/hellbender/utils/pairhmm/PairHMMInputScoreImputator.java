package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Common interface for pair-hmm score calculators.
 */
public interface PairHMMInputScoreImputator {

    /**
     * Given a read returns an score imputation that in turn can be queried for the pair-HMM match, insertion, deletion scores
     * at each position of the read.
     * @param read the target read.
     * @return never {@code null}.
     */
    PairHMMInputScoreImputation impute(final GATKRead read);

}
