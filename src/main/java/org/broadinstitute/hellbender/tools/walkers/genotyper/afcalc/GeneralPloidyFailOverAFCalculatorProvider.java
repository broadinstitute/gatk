package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Provider that defaults to the general ploidy implementation when the preferred one does not handle the required
 * ploidy.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class GeneralPloidyFailOverAFCalculatorProvider extends AFCalculatorProvider {

    private final AFCalculator preferred;
    private final AFCalculatorImplementation preferredImplementation;
    private final AFCalculator failOver;

    /**
     * Creates a new AF calculator provider given the genotyping arguments and logger reference.
     * @param genotypeArgs genotyping parameter collection.
     * @throws IllegalArgumentException if {@code genotypeArgs} is {@code null}.
     */
    public GeneralPloidyFailOverAFCalculatorProvider(final GenotypeCalculationArgumentCollection genotypeArgs) {
        Utils.nonNull(genotypeArgs);
        preferredImplementation = AFCalculatorImplementation.bestValue(genotypeArgs.samplePloidy,genotypeArgs.MAX_ALTERNATE_ALLELES, null);
        preferred = preferredImplementation.newInstance();
        failOver = AFCalculatorImplementation.EXACT_GENERAL_PLOIDY.newInstance();
    }

    /**
     * {@inheritDoc}
     * @param ploidy {@inheritDoc}
     * @param maximumAlternativeAlleles {@inheritDoc}
     * @return {@inheritDoc}
     */
    @Override
    public AFCalculator getInstance(final int ploidy, final int maximumAlternativeAlleles) {
        return preferredImplementation.usableForParams(ploidy,maximumAlternativeAlleles) ? preferred : failOver;
    }

}
