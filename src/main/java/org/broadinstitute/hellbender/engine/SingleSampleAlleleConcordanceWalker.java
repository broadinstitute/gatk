package org.broadinstitute.hellbender.engine;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;

/**
 * Base class for single sample concordance walkers which define concordance as agreement between records on the ref and first alt allele.
 *
 * Subclasses must implement the {@link #apply} method to process each variant
 * and may optionally implement {@link #onTraversalStart} and/or {@link #onTraversalSuccess}.
 */
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

    protected boolean doAltAllelesMatch(final VariantContext truth, final VariantContext eval) {
        final boolean sameRefAllele = truth.getReference().equals(eval.getReference());

        // We make sure that the truth has at least one alt allele.
        // If it does, we pick the first for comparison:
        final boolean containsAltAllele =
                (truth.getAlternateAlleles().size() == eval.getAlternateAlleles().size()) &&
                        ((truth.getAlternateAlleles().size() > 0) &&
                                eval.getAlternateAlleles().contains(truth.getAlternateAllele(0)));

        return sameRefAllele && containsAltAllele;
    }

    @Override
    protected boolean shouldVariantsBeMatched(final VariantContext truth, final VariantContext eval) {
        return doAltAllelesMatch(truth, eval);
    }
}
