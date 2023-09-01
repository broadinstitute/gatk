package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

import java.util.ArrayList;
import java.util.List;

/**
 * Stratifies by the FILTER status (PASS, FAIL) of the eval records
 */
public class Filter extends VariantStratifier {
    public Filter(VariantEvalEngine engine) {
        super(engine);

        states.add("called");
        states.add("filtered");
        states.add("raw");
    }

    @Override
    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        final ArrayList<Object> relevantStates = new ArrayList<>();

        relevantStates.add("raw");
        if (eval != null) {
            relevantStates.add(eval.isFiltered() ? "filtered" : "called");
        }

        return relevantStates;
    }
}
