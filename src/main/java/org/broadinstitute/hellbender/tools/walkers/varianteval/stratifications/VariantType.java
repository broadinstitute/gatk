package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.Collections;
import java.util.List;

/**
 * Stratifies the eval variants by their type (SNP, INDEL, ETC)
 */
public class VariantType extends VariantStratifier {
    @Override
    public void initialize() {
        for (VariantContext.Type t : VariantContext.Type.values())
            states.add(t.toString());
    }

    public List<Object> getRelevantStates(ReferenceContext referenceContext, ReadsContext readsContext, FeatureContext featureContext, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        return eval == null ? Collections.emptyList() : Collections.singletonList((Object)eval.getType().toString());
    }

}
