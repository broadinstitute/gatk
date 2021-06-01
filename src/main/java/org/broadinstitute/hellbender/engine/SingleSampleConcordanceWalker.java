package org.broadinstitute.hellbender.engine;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.validation.ConcordanceState;

/**
 * Base class for single sample concordance walkers, which process one variant at a time from one or more sources of variants,
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

    protected SingleSampleTruthVersusEval convertToSingleSampleTruthVersusEval(final TruthVersusEval truthVersusEval) {
        final Genotype truthGenotype = truthVersusEval.hasTruth() ? truthVersusEval.getTruth().getGenotype(0) : null;
        final Genotype evalGenotype = truthVersusEval.hasEval() ? truthVersusEval.getEval().getGenotype(0) : null;
        final ConcordanceState concordanceState = getConcordanceState(truthGenotype, evalGenotype, truthVersusEval.hasEval() && truthVersusEval.getEval().isFiltered());

        return new SingleSampleTruthVersusEval(truthVersusEval.getTruth(), truthVersusEval.getEval(), concordanceState);
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

        public SingleSampleTruthVersusEval(final VariantContext truth, final VariantContext eval,
                                final ConcordanceState concordanceState) {
            super(truth, eval);
            this.concordanceState = concordanceState;
        }

        public ConcordanceState getConcordance() {
            return concordanceState;
        }

    }

}
