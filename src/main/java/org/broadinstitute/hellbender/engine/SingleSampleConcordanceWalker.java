package org.broadinstitute.hellbender.engine;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;

import java.util.ArrayDeque;
import java.util.HashSet;
import java.util.List;
import java.util.Queue;
import java.util.Set;
import java.util.function.BiFunction;

/**
 * Base class for concordance walkers, which process one variant at a time from one or more sources of variants,
 * with optional contextual information from a reference, sets of reads, and/or supplementary sources of
 * Features.
 *
 * Subclasses must implement the {@link #apply} and {@link #areVariantsAtSameLocusConcordant} methods to process each variant
 * and may optionally implement {@link #onTraversalStart} and/or {@link #onTraversalSuccess}.
 *
 * Created by Takuto Sato 1/30/17, abstractified by David Benjamin on 2/22/17.
 * {@link #onTraversalStart}, {@link #onTraversalSuccess} and/or {@link #closeTool}.
 */
public abstract class SingleSampleConcordanceWalker extends AbstractConcordanceWalker {

    @Override
    protected void apply(final TruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext) {
        apply(convertToSingleSampleTruthVersusEval(truthVersusEval), readsContext, refContext);
    }

    // the primary work of the walker.  Must be overridden in implementing classes.
    protected abstract void apply(final SingleSampleTruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext);

    private SingleSampleTruthVersusEval convertToSingleSampleTruthVersusEval(final TruthVersusEval truthVersusEval) {
        if (!truthVersusEval.hasEval()) {
            return SingleSampleTruthVersusEval.falseNegative(truthVersusEval.getTruth());
        } else if (!truthVersusEval.hasTruth()) {
            final VariantContext eval = truthVersusEval.getEval();
            return eval.isFiltered() ? SingleSampleTruthVersusEval.filteredTrueNegative(eval) : SingleSampleTruthVersusEval.falsePositive(truthVersusEval.getEval());
        } else if (truthVersusEval.getEval().isFiltered()) {
            return SingleSampleTruthVersusEval.filteredFalseNegative(truthVersusEval.getTruth(), truthVersusEval.getEval());
        } else  {
            //we only match two variants together if they are concordant
            return SingleSampleTruthVersusEval.truePositive(truthVersusEval.getTruth(), truthVersusEval.getEval());
        }
    }

    // override this to customize the degree of agreement required to call a true positive.  Sometimes, for example, we may want
    // just a single alt allele to agree and sometime we may require all alt alleles to match
    protected abstract boolean areVariantsAtSameLocusConcordant(final VariantContext truth, final VariantContext eval);

    //match variant contexts based on areVariantsAtSameLocusConcordant.  This way we will only match up two variants if they are concordant
    @Override
    protected boolean shouldVariantsBeMatched(final VariantContext truth, final VariantContext eval) {
        return areVariantsAtSameLocusConcordant(truth, eval);
    }
    /**
     * store a truth vc in case of a false negative, an eval vc in case of a false positive, or a concordance pair of
     * truth and eval in case of a true positive.
     */
    protected static class SingleSampleTruthVersusEval extends TruthVersusEval {
        private final ConcordanceState concordanceState;

        private SingleSampleTruthVersusEval(final VariantContext truth, final VariantContext eval,
                                final ConcordanceState concordanceState) {
            super(truth, eval);
            this.concordanceState = concordanceState;
        }

        public static SingleSampleTruthVersusEval falseNegative(final VariantContext truth) {
            return new SingleSampleTruthVersusEval(truth, null, ConcordanceState.FALSE_NEGATIVE);
        }

        public static SingleSampleTruthVersusEval falsePositive(final VariantContext eval) {
            return new SingleSampleTruthVersusEval(null, eval, ConcordanceState.FALSE_POSITIVE);
        }

        public static SingleSampleTruthVersusEval truePositive(final VariantContext truth, final VariantContext eval) {
            return new SingleSampleTruthVersusEval(truth, eval, ConcordanceState.TRUE_POSITIVE);
        }

        public static SingleSampleTruthVersusEval filteredFalseNegative(final VariantContext truth, final VariantContext eval) {
            return new SingleSampleTruthVersusEval(truth, eval, ConcordanceState.FILTERED_FALSE_NEGATIVE);
        }

        public static SingleSampleTruthVersusEval filteredTrueNegative(final VariantContext eval) {
            return new SingleSampleTruthVersusEval(null, eval, ConcordanceState.FILTERED_TRUE_NEGATIVE);
        }

        public ConcordanceState getConcordance() {
            return concordanceState;
        }

    }

}
