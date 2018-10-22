package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Stratifies the evaluation by each contig in the reference sequence
 */
public class Contig extends VariantStratifier {
    @Override
    public void initialize() {
        states.addAll(getVariantEvalWalker().getContigNames());
        states.add("all");
    }

    @Override
    public List<Object> getRelevantStates(ReferenceContext referenceContext, ReadsContext readsContext, FeatureContext featureContext, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        if (eval != null) {
            return Arrays.asList("all", eval.getContig());
        } else {
            return Collections.emptyList();
        }
    }
}
