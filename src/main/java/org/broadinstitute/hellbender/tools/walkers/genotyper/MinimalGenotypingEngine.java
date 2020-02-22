package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
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
public final class MinimalGenotypingEngine extends GenotypingEngine<StandardCallerArgumentCollection> {

    private final DragstrParams dragstrParams;

    private ReferenceContext referenceContext;

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
        dragstrParams = null;
    }

    /**
     * Creates a new genotyping engine given the configuration parameters and the targeted set of samples
     *
     * @param configuration the configuration.
     * @param samples list of samples
     * @param doAlleleSpecificCalcs Whether to calculate genotyping annotations needed for allele specific annotations
     */
    public MinimalGenotypingEngine(final StandardCallerArgumentCollection configuration, final SampleList samples, boolean doAlleleSpecificCalcs, final DragstrParams dragstrParams) {
        super(configuration, samples, doAlleleSpecificCalcs);
        this.dragstrParams = dragstrParams;
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

