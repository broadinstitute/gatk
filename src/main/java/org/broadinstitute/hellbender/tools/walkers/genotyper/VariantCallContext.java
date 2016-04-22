package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * Useful helper class to communicate the results of calculateGenotype to framework
 */
public final class VariantCallContext extends VariantContext {

    private static final long serialVersionUID = 1L;

    // Was the site called confidently, either reference or variant?
    private boolean confidentlyCalled = false;

    // Should this site be emitted?
    private boolean shouldEmit = true;

    VariantCallContext(VariantContext vc, boolean confidentlyCalledP) {
        super(vc);
        this.confidentlyCalled = confidentlyCalledP;
    }

    VariantCallContext(VariantContext vc, boolean confidentlyCalledP, boolean shouldEmit) {
        super(vc);
        this.confidentlyCalled = confidentlyCalledP;
        this.shouldEmit = shouldEmit;
    }

    /* these methods are only implemented for GENOTYPE_GIVEN_ALLELES MODE */
    //todo -- expand these methods to all modes

    /**
     *
     * @param callConfidenceThreshold the Unified Argument Collection STANDARD_CONFIDENCE_FOR_CALLING
     * @return true if call was confidently ref
     */
    public boolean isCalledRef(double callConfidenceThreshold) {
        return (confidentlyCalled && (getPhredScaledQual() < callConfidenceThreshold));
    }

    /**
     *
     * @param callConfidenceThreshold the Unified Argument Collection STANDARD_CONFIDENCE_FOR_CALLING
     * @return true if call was confidently alt
     */
    public boolean isCalledAlt(double callConfidenceThreshold) {
        return (confidentlyCalled && (getPhredScaledQual() > callConfidenceThreshold));
    }

}