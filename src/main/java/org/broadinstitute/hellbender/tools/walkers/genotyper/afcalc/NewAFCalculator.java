package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAlleleCounts;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.hellbender.utils.Dirichlet;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class NewAFCalculator extends AFCalculator {
    private static final GenotypeLikelihoodCalculators GL_CALCS = new GenotypeLikelihoodCalculators();
    private static final double THRESHOLD_FOR_ALLELE_COUNT_CONVERGENCE = 0.1;
    private static final int HOM_REF_GENOTYPE_INDEX = 0;

    private final double refPseudocount;
    private final double snpPseudocount;
    private final double indelPseudocount;
    private final int defaultPloidy;


    public NewAFCalculator(final double refPseudocount, final double snpPseudocount, final double indelPseudocount, final int defaultPloidy) {
        this.refPseudocount = refPseudocount;
        this.snpPseudocount = snpPseudocount;
        this.indelPseudocount = indelPseudocount;
        this.defaultPloidy = defaultPloidy;
    }

    public AFCalculationResult getLog10PNonRef(final VariantContext vc) {
        return getLog10PNonRef(vc, defaultPloidy, 0, null);
    }
    //TODO: this should be a class of static methods once the old AFCalculator is gone.
    /**
     * Compute the probability of the alleles segregating given the genotype likelihoods of the samples in vc
     *
     * @param vc the VariantContext holding the alleles and sample information.  The VariantContext
     *           must have at least 1 alternative allele
     * @param refSnpIndelPseudocounts a total hack.  A length-3 vector containing Dirichlet prior pseudocounts to
     *                                be given to ref, alt SNP, and alt indel alleles.  Hack won't be necessary when we destroy the old AF calculators
     * @return result (for programming convenience)
     */
    @Override
    public AFCalculationResult getLog10PNonRef(final VariantContext vc, final int defaultPloidy, final int maximumAlternativeAlleles, final double[] refSnpIndelPseudocounts) {
        Utils.nonNull(vc, "VariantContext cannot be null");
        final int numAlleles = vc.getNAlleles();
        final List<Allele> alleles = vc.getAlleles();
        Utils.validateArg( numAlleles > 1, () -> "VariantContext has only a single reference allele, but getLog10PNonRef requires at least one at all " + vc);

        final double[] priorPseudocounts = alleles.stream()
                .mapToDouble(a -> a.isReference() ? refPseudocount : (a.length() > 1 ? snpPseudocount : indelPseudocount)).toArray();

        final Dirichlet prior = new Dirichlet(priorPseudocounts);
        double[] alleleCounts = new double[numAlleles];
        final double flatLog10AlleleFrequency = -MathUtils.log10(numAlleles); // log10(1/numAlleles)
        double[] log10AlleleFrequencies = new IndexRange(0, numAlleles).mapToDouble(n -> flatLog10AlleleFrequency);
        double alleleCountsMaximumDifference = Double.POSITIVE_INFINITY;

        while (alleleCountsMaximumDifference > THRESHOLD_FOR_ALLELE_COUNT_CONVERGENCE) {
            final double[] newAlleleCounts = effectiveAlleleCounts(vc, log10AlleleFrequencies);
            alleleCountsMaximumDifference = Arrays.stream(MathArrays.ebeSubtract(alleleCounts, newAlleleCounts)).map(Math::abs).max().getAsDouble();
            alleleCounts = newAlleleCounts;

            // first iteration uses flat prior in order to avoid local minimum where the prior + no pseudocounts gives such a low
            // effective allele frequency that it overwhelms the genotype likelihood of a real variant
            // basically, we want a chance to get non-zero pseudocounts before using a prior that's biased against a variant
            log10AlleleFrequencies = new Dirichlet(prior, alleleCounts).log10MeanWeights();
        }

        double[] log10POfZeroCountsByAllele = new double[numAlleles];
        double log10PNoVariant = 0;

        for (final Genotype g : vc.getGenotypes()) {
            if (!g.hasLikelihoods()) {
                continue;
            }
            final int ploidy = g.getPloidy() == 0 ? defaultPloidy : g.getPloidy();
            final GenotypeLikelihoodCalculator glCalc = GL_CALCS.getInstance(ploidy, numAlleles);
            final double[] genotypePosteriors = normalizedGenotypePosteriors(g, glCalc, log10AlleleFrequencies);

            //the total probability
            log10PNoVariant += Math.log10(genotypePosteriors[HOM_REF_GENOTYPE_INDEX]);

            // per allele non-log space probabilities of zero counts for this sample


            // for each allele calculate the total probability of genotypes containing at least one copy of the allele
            double[] pOfNonZeroAltAlleles = new double[numAlleles];
            new IndexRange(0, glCalc.genotypeCount()).forEach(genotypeIndex ->
                    glCalc.genotypeAlleleCountsAt(genotypeIndex).forEachAlleleIndexAndCount((alleleIndex, count) ->
                            pOfNonZeroAltAlleles[alleleIndex] += genotypePosteriors[genotypeIndex]));

            new IndexRange(0, numAlleles).forEach(a -> log10POfZeroCountsByAllele[a] += Math.log10(1 - pOfNonZeroAltAlleles[a]));
        }

        // unfortunately AFCalculationResult expects integers for the MLE.  We really should emit the EM no-integer values
        // which are valuable (eg in CombineGVCFs) as the sufficient statistics of the Dirichlet posterior on allele frequencies
        final int[] integerAlleleCounts = Arrays.stream(alleleCounts).mapToInt(x -> (int) Math.round(x)).toArray();
        final int[] integerAltAlleleCounts = Arrays.copyOfRange(integerAlleleCounts, 1, numAlleles);

        //skip the ref allele (index 0)
        final Map<Allele, Double> log10PRefByAllele = IntStream.range(1, numAlleles).boxed()
                .collect(Collectors.toMap(alleles::get, a -> log10POfZeroCountsByAllele[a]));

        // we compute posteriors here and don't have the same prior that AFCalculationResult expects.  Therefore, we
        // give it our posterior as its "likelihood" along with a flat dummy prior
        final double[] dummyFlatPrior = {-1e-10, -1e-10};   //TODO: HACK must be negative for AFCalcResult
        final double[] log10PosteriorOfNoVariantYesVariant = {log10PNoVariant, Math.log10(1 - Math.pow(10, log10PNoVariant))};

        return new AFCalculationResult(integerAltAlleleCounts, alleles, log10PosteriorOfNoVariantYesVariant, dummyFlatPrior, log10PRefByAllele);
    }

    private double[] effectiveAlleleCounts(final VariantContext vc, final double[] log10AlleleFrequencies) {
        final int numAlleles = vc.getNAlleles();
        Utils.validateArg(numAlleles == log10AlleleFrequencies.length, "number of alleles inconsistent");
        double[] result = new double[numAlleles];
        for (final Genotype g : vc.getGenotypes()) {
            if (!g.hasLikelihoods()) {
                continue;
            }
            final GenotypeLikelihoodCalculator glCalc = GL_CALCS.getInstance(g.getPloidy(), numAlleles);
            final double[] genotypePosteriors = normalizedGenotypePosteriors(g, glCalc, log10AlleleFrequencies);

            new IndexRange(0, glCalc.genotypeCount()).forEach(genotypeIndex ->
                glCalc.genotypeAlleleCountsAt(genotypeIndex).forEachAlleleIndexAndCount((alleleIndex, count) ->
                        result[alleleIndex] += genotypePosteriors[genotypeIndex] * count));
        }
        return result;
    }

    private static double[] normalizedGenotypePosteriors(final Genotype g, final GenotypeLikelihoodCalculator glCalc, final double[] log10AlleleFrequencies) {
        final double[] log10Likelihoods = g.getLikelihoods().getAsVector();
        final double[] unnormalizedLog10Likelihoods = new IndexRange(0, glCalc.genotypeCount()).mapToDouble(genotypeIndex -> {
            final GenotypeAlleleCounts gac = glCalc.genotypeAlleleCountsAt(genotypeIndex);
            return gac.log10CombinationCount() + log10Likelihoods[genotypeIndex]
                    + gac.sumOverAlleleIndicesAndCounts((index, count) -> count * log10AlleleFrequencies[index]);
        });
        return MathUtils.normalizeFromLog10(unnormalizedLog10Likelihoods);
    }

    @Override   //Note: unused
    protected AFCalculationResult getResultFromFinalState(final VariantContext vc, final double[] priors, final StateTracker st) { return null; }

    @Override//Note: unused
    protected AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy,
                                                               final double[] priors, final StateTracker st) { return null; }

    @Override   //Note: unused
    protected StateTracker getStateTracker(final boolean reset, final int maximumAlternativeAlleleCount) { return null; }

}
