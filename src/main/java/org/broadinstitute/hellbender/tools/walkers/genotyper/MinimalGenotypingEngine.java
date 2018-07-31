package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculatorProvider;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;


/**
 * A stripped-down version of the former UnifiedGenotyper's genotyping strategy implementation,
 * used only by the HaplotypeCaller for its isActive() determination. Should not be used for
 * any other purpose!
 */
public final class MinimalGenotypingEngine extends GenotypingEngine<UnifiedArgumentCollection> {

    /**
     * Creates a new unified genotyping given the UG configuration parameters and the targeted set of samples
     *
     * @param configuration the UG configuration.
     * @param samples list of samples
     */
    public MinimalGenotypingEngine(final UnifiedArgumentCollection configuration, final SampleList samples,
                                   final AFCalculatorProvider afCalculatorProvider) {
        this(configuration, samples, afCalculatorProvider, false);
    }

    /**
     * Creates a new unified genotyping given the UG configuration parameters and the targeted set of samples
     *
     * @param configuration the UG configuration.
     * @param samples list of samples
     * @param doAlleleSpecificCalcs Whether to calculate genotyping annotations needed for allele specific annotations
     */
    public MinimalGenotypingEngine(final UnifiedArgumentCollection configuration, final SampleList samples,
                                    final AFCalculatorProvider afCalculatorProvider, boolean doAlleleSpecificCalcs ) {
        super(configuration, samples, afCalculatorProvider, doAlleleSpecificCalcs);

        if ( configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES ) {
            throw new UserException("GENOTYPE_GIVEN_ALLELES mode not supported in the MinimalGenotypingEngine");
        } else if ( configuration.GLmodel != GenotypeLikelihoodsCalculationModel.SNP ) {
            throw new UserException("Only the diploid SNP model is supported in the MinimalGenotypingEngine");
        } else if ( configuration.COMPUTE_SLOD ) {
            throw new UserException("--computeSLOD not supported in the MinimalGenotypingEngine");
        }
    }

    @Override
    protected boolean forceKeepAllele(final Allele allele) {
        return configuration.genotypingOutputMode == GenotypingOutputMode.GENOTYPE_GIVEN_ALLELES || configuration.annotateAllSitesWithPLs;
    }

    @Override
    protected String callSourceString() {
        return "UG_call";
    }

    @Override
    protected boolean forceSiteEmission() {
        return configuration.outputMode == OutputMode.EMIT_ALL_SITES;
    }
}

