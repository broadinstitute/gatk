package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Stratifies by whether a site in in the list of known RODs (e.g., dbsnp by default)
 */
public class Novelty extends VariantStratifier implements StandardStratification {
    // needs the variant contexts and known names
    private List<FeatureInput<VariantContext>> knowns;

    private final static List<Object> KNOWN_STATES = Arrays.asList((Object)"all", (Object)"known");
    private final static List<Object> NOVEL_STATES = Arrays.asList((Object)"all", (Object)"novel");

    public Novelty(VariantEvalEngine engine) {
        super(engine);

        states.addAll(Arrays.asList("all", "known", "novel"));
        knowns = getEngine().getKnowns();
    }

    @Override
    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        if (eval != null) {
            //NOTE: this is limiting to matching start position, comparable to GATK3.  Unsure if we shoud carry behavior that forward
            final Collection<VariantContext> knownComps = knowns.stream().map(context::getVariantsForFeature).flatMap(List::stream).collect(Collectors.toCollection(HashSet::new));
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