package org.broadinstitute.hellbender.tools.haplotypecaller;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;

import java.util.List;

public interface GenotypePriorCalculator {
    double[] getLog10Priors(final GenotypeLikelihoodCalculator lkCalculator, final List<Allele> alleles);

}
