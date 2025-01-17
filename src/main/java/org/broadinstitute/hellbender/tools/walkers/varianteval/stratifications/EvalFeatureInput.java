package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalArgumentCollection;
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

        VariantEvalArgumentCollection args = engine.getVariantEvalArgs();
        for ( FeatureInput<VariantContext> fi : args.getEvals() ) {
            states.add(engine.getNameForInput(fi));
            if ( args.isMergeEvals() )
                break;
        }
    }

    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        return Collections.singletonList(evalName);
    }
}
