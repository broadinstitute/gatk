package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications.manager.Stratifier;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

public abstract class VariantStratifier implements Comparable<VariantStratifier>, Stratifier<Object> {
    final private String name;
    final private VariantEvalEngine engine;

    final protected ArrayList<Object> states = new ArrayList<>();

    public VariantStratifier(VariantEvalEngine engine) {
        this.name = this.getClass().getSimpleName();
        this.engine = engine;
    }

    public VariantEvalEngine getEngine() {
        return engine;
    }

    // Subclasses can override to validate arguments
    public void validateArgs() {

    }

    public abstract List<Object> getRelevantStates(VariantEvalContext context, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String familyName);

    // -------------------------------------------------------------------------------------
    //
    // final capabilities
    //
    // -------------------------------------------------------------------------------------

    public final int compareTo(VariantStratifier o1) {
        return this.getName().compareTo(o1.getName());
    }

    @Override
    public String toString() {
        return getName();
    }

    public final String getName() {
        return name;
    }
    
    public String getFormat() { return "%s"; }
    
    public final ArrayList<Object> getAllStates() {
        return states;
    }


    /**
     * The way for a stratifier to specify that it's incompatible with specific evaluations.  For
     * example, VariantSummary includes a per-sample metric, and so cannot be used safely with Sample
     * or AlleleCount stratifications as this introduces an O(n^2) memory and cpu cost.
     *
     * @return the set of VariantEvaluators that cannot be active with this Stratification
     */
    public Set<Class<? extends VariantEvaluator>> getIncompatibleEvaluators() {
        return Collections.emptySet();
    }
}
