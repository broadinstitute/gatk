package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantSummary;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

import java.util.*;

/**
 * Stratifies the eval RODs by each sample in the eval ROD.
 *
 * This allows the system to analyze each sample separately.  Since many evaluations
 * only consider non-reference sites, stratifying by sample results in meaningful
 * calculations for CompOverlap
 */
public class Sample extends VariantStratifier {
    public Sample(VariantEvalEngine engine) {
        super(engine);

        states.addAll(getEngine().getSampleNamesForStratification());
    }

    @Override
    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        return Collections.singletonList(sampleName);
    }

    @Override
    public Set<Class<? extends VariantEvaluator>> getIncompatibleEvaluators() {
        return new HashSet<>(Arrays.asList(VariantSummary.class));
    }
}
