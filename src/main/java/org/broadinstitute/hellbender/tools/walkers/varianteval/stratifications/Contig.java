package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Stratifies the evaluation by each contig in the reference sequence. Note: if the user supplies custom intervals, it will defer to these rather than the full sequence dictionary
 */
public class Contig extends VariantStratifier {
    public Contig(VariantEvalEngine engine) {
        super(engine);

        states.addAll(getEngine().getContigNames());
        states.add("all");
    }

    @Override
    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        if (eval != null) {
            return Arrays.asList("all", eval.getContig());
        } else {
            return Collections.emptyList();
        }
    }
}
