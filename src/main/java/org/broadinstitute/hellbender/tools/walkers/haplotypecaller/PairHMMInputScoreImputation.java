package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

public interface PairHMMInputScoreImputation {

    byte[] delOpenPenalties();

    byte[] insOpenPenalties();

    byte[] gapContinuationPenalties();

}
