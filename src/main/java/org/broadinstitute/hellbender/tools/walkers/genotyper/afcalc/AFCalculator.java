package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;


/**
 * Generic interface for calculating the probability of alleles segregating given priors and genotype likelihoods
 */
public abstract class AFCalculator {

    protected static final Logger logger = LogManager.getLogger(AFCalculator.class);

    private StateTracker stateTracker;

    /**
     * Compute the probability of the alleles segregating given the genotype likelihoods of the samples in vc
     *
     * @param vc the VariantContext holding the alleles and sample information.  The VariantContext
     *           must have at least 1 alternative allele
     * @param log10AlleleFrequencyPriors a prior vector nSamples x 2 in length indicating the Pr(AF = i)
     * @return result (for programming convenience)
     */
    public AFCalculationResult getLog10PNonRef(final VariantContext vc, final int defaultPloidy, final int maximumAlternativeAlleles, final double[] log10AlleleFrequencyPriors) {
        Utils.nonNull(vc, "VariantContext cannot be null");
        Utils.nonNull(log10AlleleFrequencyPriors, "priors vector cannot be null");
        if ( vc.getNAlleles() == 1 ) {
            throw new IllegalArgumentException("VariantContext has only a single reference allele, but getLog10PNonRef requires at least one at all " + vc);
        }

        // reset the result, so we can store our new result there
        final StateTracker stateTracker = getStateTracker(true,maximumAlternativeAlleles);

        final VariantContext vcWorking = reduceScope(vc,defaultPloidy, maximumAlternativeAlleles);

        return computeLog10PNonRef(vcWorking, defaultPloidy, log10AlleleFrequencyPriors, stateTracker);
    }

    /**
     * Convert the final state of the state tracker into our result as an AFCalculationResult
     *
     * Assumes that stateTracker has been updated accordingly
     *
     * @param vcWorking the VariantContext we actually used as input to the calc model (after reduction)
     * @param log10AlleleFrequencyPriors the priors by AC vector
     * @return a AFCalculationResult describing the result of this calculation
     */
    protected AFCalculationResult getResultFromFinalState(final VariantContext vcWorking, final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker) {
        Utils.nonNull(vcWorking, "vcWorking cannot be null");
        Utils.nonNull(log10AlleleFrequencyPriors, "log10AlleleFrequencyPriors cannot be null");

        stateTracker.setAllelesUsedInGenotyping(vcWorking.getAlleles());
        return stateTracker.toAFCalculationResult(log10AlleleFrequencyPriors);
    }

    // ---------------------------------------------------------------------------
    //
    // Abstract methods that should be implemented by concrete implementations
    // to actually calculate the AF
    //
    // ---------------------------------------------------------------------------

    /**
     * Look at VC and perhaps return a new one of reduced complexity, if that's necessary
     *
     * Used before the call to computeLog10PNonRef to simply the calculation job at hand,
     * if vc exceeds bounds.  For example, if VC has 100 alt alleles this function
     * may decide to only genotype the best 2 of them.
     *
     * @param vc the initial VC provided by the caller to this AFcalculation
     * @return a potentially simpler VC that's more tractable to genotype
     */
    protected abstract VariantContext reduceScope(final VariantContext vc, final int defaultPloidy, final int maximumAlternativeAlleles);

    /**
     * Actually carry out the log10PNonRef calculation on vc, storing results in results
     *
     * @param vc                                variant context with alleles and genotype likelihoods,
     *                                          must have at least one alt allele
     * @param log10AlleleFrequencyPriors        priors
     * @return a AFCalcResult object describing the results of this calculation
     */
    protected abstract AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy,
                                                        final double[] log10AlleleFrequencyPriors, final StateTracker stateTracker);

    /**
     * Subset VC to the just allelesToUse, updating genotype likelihoods
     *
     * Must be overridden by concrete subclasses
     *
     * @param vc                                variant context with alleles and genotype likelihoods
     * @param defaultPloidy                     default ploidy to assume in case {@code vc} does not indicate it for a sample.
     * @param allelesToUse                      alleles to subset
     * @param assignGenotypes
     * @return GenotypesContext object
     */
    public abstract GenotypesContext subsetAlleles(final VariantContext vc,
                                                   final int defaultPloidy,
                                                   final List<Allele> allelesToUse,
                                                   final boolean assignGenotypes);

    /**
     * Retrieves the state tracker.
     *
     * <p>
     *     The tracker will be reset if so requested or if it needs to be resized due to an increase in the
     *     maximum number of alleles is must be able to handle.
     * </p>
     *
     * @param reset make sure the tracker is reset.
     * @param maximumAlternativeAlleleCount the maximum alternative allele count it must be able to handle. Has no effect if
     *                                     the current tracker is able to handle that number.
     *
     * @return never {@code null}
     */
    protected StateTracker getStateTracker(final boolean reset, final int maximumAlternativeAlleleCount) {
        if (stateTracker == null) {
            stateTracker = new StateTracker(maximumAlternativeAlleleCount);
        } else if (reset) {
            stateTracker.reset(maximumAlternativeAlleleCount);
        } else {
            stateTracker.ensureMaximumAlleleCapacity(maximumAlternativeAlleleCount);
        }
        return stateTracker;
    }

    /**
     * Please don't use this method in production.
     */
    @VisibleForTesting
    int getAltAlleleCountOfMAP(final int allele) {
        return getStateTracker(false,allele + 1).getAlleleCountsOfMAP()[allele];
    }

}