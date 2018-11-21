package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * Stratifies by whether a site in in the list of known RODs (e.g., dbsnp by default)
 */
public class Novelty extends VariantStratifier implements StandardStratification {
    // needs the variant contexts and known names
    private List<FeatureInput<VariantContext>> knowns;

    private final static List<Object> KNOWN_STATES = Arrays.asList((Object)"all", (Object)"known");
    private final static List<Object> NOVEL_STATES = Arrays.asList((Object)"all", (Object)"novel");

    @Override
    public void initialize() {
        states.addAll(Arrays.asList("all", "known", "novel"));
        knowns = getVariantEvalWalker().getKnowns();
    }

    public List<Object> getRelevantStates(ReferenceContext referenceContext, ReadsContext readsContext, FeatureContext featureContext, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        if (eval != null) {
            //NOTE: this is limiting to matching start position, comparable to GATK3.  Unsure if we shoud carry behavior that forward
            final Collection<VariantContext> knownComps = featureContext.getValues(knowns, eval.getStart());
            for ( final VariantContext c : knownComps ) {
                // loop over sites, looking for something that matches the type eval
                if ( eval.getType() == c.getType() || eval.getType() == VariantContext.Type.NO_VARIATION ) {
                    return KNOWN_STATES;
                }
            }
        } 
        
        return NOVEL_STATES;
    }
}