package org.broadinstitute.hellbender.utils.pairhmm;

public interface PairHMMInputScoreImputation {

    byte[] delOpenPenalties();

    byte[] insOpenPenalties();

    byte[] gapContinuationPenalties();

}
