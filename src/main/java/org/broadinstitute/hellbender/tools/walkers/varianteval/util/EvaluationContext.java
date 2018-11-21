package org.broadinstitute.hellbender.tools.walkers.varianteval.util;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEval;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.manager.StratificationManager;

import java.util.ArrayList;
import java.util.Set;
import java.util.TreeSet;

public final class EvaluationContext {
    // NOTE: must be hashset to avoid O(log n) cost of iteration in the very frequently called apply function
    final VariantEval walker;
    private final ArrayList<VariantEvaluator> evaluationInstances;
    private final Set<Class<? extends VariantEvaluator>> evaluationClasses;

    public EvaluationContext(final VariantEval walker, final Set<Class<? extends VariantEvaluator>> evaluationClasses) {
        this(walker, evaluationClasses, true);
    }

    private EvaluationContext(final VariantEval walker, final Set<Class<? extends VariantEvaluator>> evaluationClasses, final boolean doInitialize) {
        this.walker = walker;
        this.evaluationClasses = evaluationClasses;
        this.evaluationInstances = new ArrayList<VariantEvaluator>(evaluationClasses.size());

        for ( final Class<? extends VariantEvaluator> c : evaluationClasses ) {
            try {
                final VariantEvaluator eval = c.newInstance();
                if ( doInitialize ) eval.initialize(walker);
                evaluationInstances.add(eval);
            } catch (InstantiationException e) {
                throw new GATKException("Unable to instantiate eval module '" + c.getSimpleName() + "'", e);
            } catch (IllegalAccessException e) {
                throw new GATKException("Illegal access error when trying to instantiate eval module '" + c.getSimpleName() + "'", e);
            }
        }
    }

    /**
     * Returns a sorted set of VariantEvaluators
     *
     * @return
     */
    public final TreeSet<VariantEvaluator> getVariantEvaluators() {
        return new TreeSet<VariantEvaluator>(evaluationInstances);
    }

    public final void apply(ReferenceContext referenceContext, ReadsContext readsContext, FeatureContext featureContext, VariantContext comp, VariantContext eval) {
        for ( final VariantEvaluator evaluation : evaluationInstances ) {
            // now call the single or paired update function
            switch ( evaluation.getComparisonOrder() ) {
                case 1:
                    if (eval != null) {
                        evaluation.update1(eval, referenceContext, readsContext, featureContext);
                    }
                    break;
                case 2:
                    evaluation.update2(eval, comp, referenceContext, readsContext, featureContext);
                    break;
                default:
                    throw new GATKException("BUG: Unexpected evaluation order " + evaluation);
            }
        }
    }

    public void combine(final EvaluationContext rhs) {
        for ( int i = 0; i < evaluationInstances.size(); i++ )
            evaluationInstances.get(i).combine(rhs.evaluationInstances.get(i));
    }

    public final static EvaluationContextCombiner COMBINER = new EvaluationContext.EvaluationContextCombiner();
    private static class EvaluationContextCombiner implements StratificationManager.Combiner<EvaluationContext> {
        @Override
        public EvaluationContext combine(EvaluationContext lhs, final EvaluationContext rhs) {
            if ( lhs == null )
                lhs = new EvaluationContext(rhs.walker, rhs.evaluationClasses, false);
            lhs.combine(rhs);
            return lhs;
        }
    }
}
