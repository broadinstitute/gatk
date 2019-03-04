package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

import java.util.List;
import java.util.Set;

/**
 * Carries the result of a call to #assignGenotypeLikelihoods
 */
public class CalledHaplotypes {
    private final List<VariantContext> calls;
    private final Set<Haplotype> calledHaplotypes;

    public CalledHaplotypes(final List<VariantContext> calls, final Set<Haplotype> calledHaplotypes) {
        this.calls = Utils.nonNull(calls, "calls cannot be null");
        this.calledHaplotypes = Utils.nonNull(calledHaplotypes, "calledHaplotypes cannot be null");
        Utils.validateArg(calls.isEmpty() == calledHaplotypes.isEmpty(), "Calls and calledHaplotypes should both be empty or both not but got calls=" + calls + " calledHaplotypes=" + calledHaplotypes);
    }

    /**
     * Get the list of calls made at this location
     * @return a non-null (but potentially empty) list of calls
     */
    public List<VariantContext> getCalls() {
        return calls;
    }

    /**
     * Get the set of haplotypes that we actually called (i.e., underlying one of the VCs in getCalls().
     * @return a non-null set of haplotypes
     */
    public Set<Haplotype> getCalledHaplotypes() {
        return calledHaplotypes;
    }
}
