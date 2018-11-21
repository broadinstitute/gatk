package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.*;

import java.util.Collections;
import java.util.List;


/**
 * Required stratification grouping output by each comp
 */
public class CompFeatureInput extends VariantStratifier implements RequiredStratification {
    @Override
    public void initialize() {
        for ( FeatureInput<VariantContext> rod : getVariantEvalWalker().getComps() ) {
            states.add(getVariantEvalWalker().getNameForInput(rod));
        }

        if (states.isEmpty()) {
            states.add("none");
        }
    }

    public List<Object> getRelevantStates(ReferenceContext referenceContext, ReadsContext readsContext, FeatureContext featureContext, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        return Collections.singletonList(compName == null ? "none" : compName);
    }
}
