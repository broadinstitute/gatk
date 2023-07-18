package org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators;

import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;

/**
 * A test evaluator which calls {@link VariantEvalEngine#getnProcessedLoci} but doesn't override {@link VariantEvaluator#requiresTerritoryToBeSpecified()}
 */
public class TestEvaluatorWhichRequiresReferenceButDoesntSayItDoes extends VariantEvaluator {
    public TestEvaluatorWhichRequiresReferenceButDoesntSayItDoes(VariantEvalEngine engine) {
        super(engine);
    }

    @Override
    public int getComparisonOrder() {
        return 1;
    }

    @Override
    public void finalizeEvaluation() {
        getEngine().getnProcessedLoci();
    }
}
