package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.utils.read.GATKRead;

public interface PairHMMInputScoreImputator {

    PairHMMInputScoreImputation impute(final GATKRead read);

}
