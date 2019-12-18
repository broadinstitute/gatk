package org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators;

import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEval;

/**
 * A test evaluator which calls {@link VariantEval#getnProcessedLoci} but doesn't override {@link VariantEvaluator#requiresTerritoryToBeSpecified()}
 */
public class TestEvaluatorWhichRequiresReferenceButDoesntSayItDoes extends VariantEvaluator {
    @Override
    public int getComparisonOrder() {
        return 1;
    }

    @Override
    public void finalizeEvaluation() {
        getWalker().getnProcessedLoci();
    }
}
