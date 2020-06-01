package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;

public interface AlleleFrequencyCalculator {

    int getPloidy();

    double[] getPriorFrequencies(AlleleList<Allele> alleleList);

    AFCalculationResult calculate(VariantContext vc);

    AFCalculationResult calculate(VariantContext vc, int defaultPloidy);

    double calculateSingleSampleBiallelicNonRefPosterior(double[] log10GenotypeLikelihoods, boolean returnZeroIfRefIsMax);
}
