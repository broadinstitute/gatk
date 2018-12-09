package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.*;

import java.util.Arrays;
import java.util.List;

/**
 * Required stratification grouping output by each eval
 */
public class EvalFeatureInput extends VariantStratifier implements RequiredStratification {
    @Override
    public void initialize() {
        for ( FeatureInput<VariantContext> rod : getVariantEvalSourceProvider().getEvals() ) {
            states.add(getVariantEvalSourceProvider().getNameForInput(rod));
            if ( getVariantEvalSourceProvider().getMergeEvals() )
                break;
        }
    }

    public List<Object> getRelevantStates(ReferenceContext referenceContext, ReadsContext readsContext, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        return Arrays.asList(evalName);
    }
}
