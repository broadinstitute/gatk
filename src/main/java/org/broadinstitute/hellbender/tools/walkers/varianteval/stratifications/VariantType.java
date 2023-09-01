package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

import java.util.Collections;
import java.util.List;

/**
 * Stratifies the eval variants by their type (SNP, INDEL, ETC)
 */
public class VariantType extends VariantStratifier {
    public VariantType(VariantEvalEngine engine) {
        super(engine);

        for (VariantContext.Type t : VariantContext.Type.values())
            states.add(t.toString());
    }

    @Override
    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        return eval == null ? Collections.emptyList() : Collections.singletonList((Object)eval.getType().toString());
    }

}
