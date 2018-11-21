package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantSummary;

import java.util.*;

/**
 * Stratifies the eval RODs by each sample in the eval ROD.
 *
 * This allows the system to analyze each sample separately.  Since many evaluations
 * only consider non-reference sites, stratifying by sample results in meaningful
 * calculations for CompOverlap
 */
public class Sample extends VariantStratifier {
    @Override
    public void initialize() {
        states.addAll(getVariantEvalWalker().getSampleNamesForStratification());
    }

    public List<Object> getRelevantStates(ReferenceContext referenceContext, ReadsContext readsContext, FeatureContext featureContext, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        return Collections.singletonList(sampleName);
    }

    @Override
    public Set<Class<? extends VariantEvaluator>> getIncompatibleEvaluators() {
        return new HashSet<>(Arrays.asList(VariantSummary.class));
    }
}
