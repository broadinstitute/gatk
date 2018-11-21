package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantSummary;

import java.util.*;

/**
 * Stratifies the eval RODs by each family in the eval ROD, as described by the pedigree.
 *
 * This allows the system to analyze each family separately.  This is particularly useful for the MendelianViolationEvaluator module.
 */
public class Family extends VariantStratifier {
    @Override
    public void initialize() {
        states.addAll(getVariantEvalWalker().getFamilyNamesForStratification());
    }

    public List<Object> getRelevantStates(ReferenceContext referenceContext, ReadsContext readsContext, FeatureContext featureContext, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String familyName) {
        return Collections.singletonList(familyName);
    }

    @Override
    public Set<Class<? extends VariantEvaluator>> getIncompatibleEvaluators() {
        return new HashSet<>(Arrays.asList(VariantSummary.class));
    }
}
