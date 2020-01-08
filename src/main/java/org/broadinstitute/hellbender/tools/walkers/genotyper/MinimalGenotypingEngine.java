package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;


/**
 * A stripped-down version of the former UnifiedGenotyper's genotyping strategy implementation,
 * used only by the HaplotypeCaller for its isActive() determination. Should not be used for
 * any other purpose!
 */
public final class MinimalGenotypingEngine extends GenotypingEngine<StandardCallerArgumentCollection> {

    /**
     * Creates a new genotyping engine given the configuration parameters and the targeted set of samples
     *
     * @param configuration the configuration.
     * @param samples list of samples
     */
    public MinimalGenotypingEngine(final StandardCallerArgumentCollection configuration, final SampleList samples) {
        this(configuration, samples, false);
    }

    /**
     * Creates a new genotyping engine given the configuration parameters and the targeted set of samples
     *
     * @param configuration the configuration.
     * @param samples list of samples
     * @param doAlleleSpecificCalcs Whether to calculate genotyping annotations needed for allele specific annotations
     */
    public MinimalGenotypingEngine(final StandardCallerArgumentCollection configuration, final SampleList samples, boolean doAlleleSpecificCalcs ) {
        super(configuration, samples, doAlleleSpecificCalcs);
    }

    @Override
    protected boolean forceKeepAllele(final Allele allele) {
        return configuration.annotateAllSitesWithPLs;
    }

    @Override
    protected String callSourceString() {
        return "UG_call";
    }
}

