package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;

public abstract class SingleSampleAlleleConcordanceWalker extends SingleSampleConcordanceWalker {

    @Override
    protected SingleSampleTruthVersusEval convertToSingleSampleTruthVersusEval(final TruthVersusEval truthVersusEval) {
        if (!truthVersusEval.hasTruth()) {
            if (truthVersusEval.getEval().isFiltered()) {
                return new SingleSampleTruthVersusEval(truthVersusEval.getTruth(), truthVersusEval.getEval(), ConcordanceState.FILTERED_TRUE_NEGATIVE);
            } else {
                return new SingleSampleTruthVersusEval(truthVersusEval.getTruth(), truthVersusEval.getEval(), ConcordanceState.FALSE_POSITIVE);
            }
        }

        if (!truthVersusEval.hasEval()) {
            return new SingleSampleTruthVersusEval(truthVersusEval.getTruth(), truthVersusEval.getEval(), ConcordanceState.FALSE_NEGATIVE);
        }

        if (truthVersusEval.getEval().isFiltered()) {
            return new SingleSampleTruthVersusEval(truthVersusEval.getTruth(), truthVersusEval.getEval(), ConcordanceState.FILTERED_FALSE_NEGATIVE);
        } else {
            return new SingleSampleTruthVersusEval(truthVersusEval.getTruth(), truthVersusEval.getEval(), ConcordanceState.TRUE_POSITIVE);
        }
    }
}
