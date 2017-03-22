package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Instantiates Exact AF calculators given the required ploidy specs.
 *
 * <p>This class might return the same instance several times
 * and so the client code might need to make sure that there are no collisions or race conditions.</p>
 */
public abstract class AFCalculatorProvider {

    /**
     * Returns a AF calculator capable to handle a particular variant-context.
     * @param variantContext the target context build.
     * @param defaultPloidy the assumed ploidy in case that there is no a GT call present to determine it.
     * @return never {@code null}
     */
    public AFCalculator getInstance(final VariantContext variantContext, final int defaultPloidy, final int maximumAltAlleles) {
        Utils.nonNull(variantContext, "variant context cannot be null");

        final int sampleCount = variantContext.getNSamples();
        if  (sampleCount == 0) {
            return getInstance(defaultPloidy, maximumAltAlleles);
        }

        final GenotypesContext genotypes = variantContext.getGenotypes();

        final Genotype firstGenotype = genotypes.get(0);
        int ploidy = firstGenotype.getPloidy();
        if (ploidy <= 0) {
            ploidy = defaultPloidy;
        }
        for (int i = 1 ; i < sampleCount; i++) {
            final Genotype genotype = genotypes.get(i);
            final int declaredPloidy = genotype.getPloidy();
            final int actualPloidy = declaredPloidy <= 0 ? defaultPloidy : declaredPloidy;
            if (actualPloidy != ploidy) {
                ploidy = AFCalculatorImplementation.UNBOUND_PLOIDY;
                break;
            }
        }
        return getInstance(ploidy, Math.min(variantContext.getNAlleles() - 1, maximumAltAlleles));
    }

    /**
     * Returns a AF calculator given the required homogeneous ploidy and maximum alt allele count.
     * @param ploidy the required ploidy.
     * @param maximumAltAlleles the maximum alt allele count.
     * @return never {@code null}
     */
    public abstract AFCalculator getInstance(final int ploidy, final int maximumAltAlleles);

}
