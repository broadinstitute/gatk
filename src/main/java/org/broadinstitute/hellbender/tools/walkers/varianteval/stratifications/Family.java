package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantSummary;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

import java.util.*;

/**
 * Stratifies the eval RODs by each family in the eval ROD, as described by the pedigree.
 *
 * This allows the system to analyze each family separately.  This is particularly useful for the MendelianViolationEvaluator module.
 */
public class Family extends VariantStratifier {
    public Family(VariantEvalEngine engine) {
        super(engine);

        states.addAll(getEngine().getFamilyNamesForStratification());
    }

    @Override
    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        return Collections.singletonList(familyName);
    }

    @Override
    public Set<Class<? extends VariantEvaluator>> getIncompatibleEvaluators() {
        return new HashSet<>(Arrays.asList(VariantSummary.class));
    }
}
