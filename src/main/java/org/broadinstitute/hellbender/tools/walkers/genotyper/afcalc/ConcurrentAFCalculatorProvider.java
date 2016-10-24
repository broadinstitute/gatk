package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * Produces independent AF calculators per thread.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public abstract class ConcurrentAFCalculatorProvider extends AFCalculatorProvider {

    private final ThreadLocal<AFCalculatorProvider> threadLocal;

    /**
     * Create a new concurrent af-calculator provider instance.
     */
    public ConcurrentAFCalculatorProvider() {
        threadLocal = new ThreadLocal<AFCalculatorProvider>() {
            @Override
            public AFCalculatorProvider initialValue() {
                return createProvider();
            }
        };
    }

    @Override
    public AFCalculator getInstance(final VariantContext vc, final int defaultPloidy, final int maxAltAlleleCount) {
        return threadLocal.get().getInstance(vc,defaultPloidy,maxAltAlleleCount);
    }


    @Override
    public AFCalculator getInstance(final int ploidy, final int maxAltAlleleCount) {
        return threadLocal.get().getInstance(ploidy, maxAltAlleleCount);
    }

    protected abstract AFCalculatorProvider createProvider();
}

