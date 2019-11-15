package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AlleleFrequencyCalculator;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.pairhmm.DragstrParams;
import org.broadinstitute.hellbender.utils.pairhmm.DragstrUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.Collections;


/**
 * A stripped-down version of the former UnifiedGenotyper's genotyping strategy implementation,
 * used only by the HaplotypeCaller for its isActive() determination. Should not be used for
 * any other purpose!
 */
public final class MinimalGenotypingEngine extends GenotypingEngine<UnifiedArgumentCollection> {

    private final DragstrParams dragstrParams;

    private ReferenceContext referenceContext;

    /**
     * Creates a new unified genotyping given the UG configuration parameters and the targeted set of samples
     *
     * @param configuration the UG configuration.
     * @param samples list of samples
     */
    public MinimalGenotypingEngine(final UnifiedArgumentCollection configuration, final SampleList samples) {
        this(configuration, samples, false, null);
    }

    /**
     * Creates a new unified genotyping given the UG configuration parameters and the targeted set of samples
     *
     * @param configuration the UG configurat
     * @param samples list of samples
     * @param doAlleleSpecificCalcs Whether to calculate genotyping annotations needed for allele specific annotations
     */
    public MinimalGenotypingEngine(final UnifiedArgumentCollection configuration, final SampleList samples, boolean doAlleleSpecificCalcs,
                                   final DragstrParams dragstrParms) {
        super(configuration, samples, doAlleleSpecificCalcs);
        this.dragstrParams = dragstrParms;
        this.referenceContext = null;

       if ( configuration.GLmodel != GenotypeLikelihoodsCalculationModel.SNP ) {
            throw new UserException("Only the diploid SNP model is supported in the MinimalGenotypingEngine");
        } else if ( configuration.COMPUTE_SLOD ) {
            throw new UserException("--computeSLOD not supported in the MinimalGenotypingEngine");
        }
    }

    @Override
    protected boolean forceKeepAllele(final Allele allele) {
        return configuration.annotateAllSitesWithPLs;
    }

    @Override
    protected String callSourceString() {
        return "UG_call";
    }

    @Override
    public VariantContext calculateGenotypes(final VariantContext vc) {
        if (dragstrParams == null || getConfiguration().genotypeArgs.dontUseDragstrPriors || !GATKVariantContextUtils.containsInlineIndel(vc)) {
            return calculateGenotypes(vc, null, Collections.emptyList());
        } else if (referenceContext != null) {
            final SimpleInterval interval = new SimpleInterval(vc.getContig(), Math.max(1, vc.getStart() - dragstrParams.maximumPeriod() * dragstrParams.maximumRepeats()), vc.getStart() - dragstrParams.maximumPeriod() * dragstrParams.maximumRepeats());

            final byte[] bases = referenceContext.getBases(interval);
            final int startOffset = vc.getStart() - interval.getStart();
            final DragstrUtils.STRSequenceAnalyzer analyzer = DragstrUtils.repeatPeriodAndCounts(bases, dragstrParams.maximumPeriod(), startOffset, startOffset + 1);
            final int period = analyzer.mostRepeatedPeriod(startOffset);
            final int repeats = analyzer.numberOfMostRepeats(startOffset);

            final AlleleFrequencyCalculator afc = dragstrParams.getAFCalculator(period, repeats, vc.getMaxPloidy(2), getConfiguration().genotypeArgs.snpHeterozygosity, getConfiguration().genotypeArgs.dragstrPriorScale);
            return  calculateGenotypes(vc, afc, Collections.emptyList());
        } else {
            return  calculateGenotypes(vc, null, Collections.emptyList());
        }
    }

    public void setReferenceContext(final ReferenceContext ref) {
        referenceContext = ref;
    }

}

