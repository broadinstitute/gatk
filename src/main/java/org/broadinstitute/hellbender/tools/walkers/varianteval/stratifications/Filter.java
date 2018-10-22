package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.ArrayList;
import java.util.List;

/**
 * Stratifies by the FILTER status (PASS, FAIL) of the eval records
 */
public class Filter extends VariantStratifier {
    @Override
    public void initialize() {
        states.add("called");
        states.add("filtered");
        states.add("raw");
    }

    public List<Object> getRelevantStates(ReferenceContext referenceContext, ReadsContext readsContext, FeatureContext featureContext, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        ArrayList<Object> relevantStates = new ArrayList<Object>();

        relevantStates.add("raw");
        if (eval != null) {
            relevantStates.add(eval.isFiltered() ? "filtered" : "called");
        }

        return relevantStates;
    }
}
