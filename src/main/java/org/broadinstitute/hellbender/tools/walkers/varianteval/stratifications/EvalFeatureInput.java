package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

import java.util.Collections;
import java.util.List;

/**
 * Required stratification grouping output by each eval
 */
public class EvalFeatureInput extends VariantStratifier implements RequiredStratification {
    public EvalFeatureInput(VariantEvalEngine engine) {
        super(engine);

        for ( FeatureInput<VariantContext> fi : getEngine().getVariantEvalArgs().getEvals() ) {
            states.add(getEngine().getNameForInput(fi));
            if ( getEngine().getVariantEvalArgs().isMergeEvals() )
                break;
        }
    }

    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        return Collections.singletonList(evalName);
    }
}
