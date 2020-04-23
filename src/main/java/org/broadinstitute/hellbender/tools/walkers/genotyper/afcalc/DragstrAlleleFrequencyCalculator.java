package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.pairhmm.DragstrParams;

public class DragstrAlleleFrequencyCalculator implements AlleleFrequencyCalculator {

    private static final double HOM_REF_PRIOR = 0;
    private static final double SNP_SIMPLE_HET_PRIOR = 34.77;
    private static final double SNP_COMPOSITE_HET_PRIOR = 69.54;
    private static final double SNP_HOMVAR_PRIOR = 37.77;

    private final double api;
    private final int defaultPloidy;

    private DragstrAlleleFrequencyCalculator(final double api, final int defaultPloidy) {
        this.api = api;
        this.defaultPloidy = defaultPloidy;
    }

    public static DragstrAlleleFrequencyCalculator makeCalculator(final DragstrParams params, final int period, final int repeats, final int defaultPloidy) {
        return new DragstrAlleleFrequencyCalculator(params.api(period, repeats), defaultPloidy);
    }

    @Override
    public int getPloidy() {
        return defaultPloidy;
    }

    @Override
    public double[] getPriorFrequencies(AlleleList<Allele> alleleList) {
        return new double[0];
    }

    @Override
    public AFCalculationResult calculate(VariantContext vc) {
        return null;
    }

    @Override
    public AFCalculationResult calculate(VariantContext vc, int defaultPloidy) {
        return null;
    }

    @Override
    public double calculateSingleSampleBiallelicNonRefPosterior(double[] log10GenotypeLikelihoods, boolean returnZeroIfRefIsMax) {
        return 0;
    }
}
