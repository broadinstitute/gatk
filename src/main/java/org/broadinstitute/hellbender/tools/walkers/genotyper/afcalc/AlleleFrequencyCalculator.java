package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAlleleCounts;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.hellbender.utils.Dirichlet;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFrequencyCalculator extends AFCalculator {
    private static final GenotypeLikelihoodCalculators GL_CALCS = new GenotypeLikelihoodCalculators();
    private static final double THRESHOLD_FOR_ALLELE_COUNT_CONVERGENCE = 0.1;
    private static final int HOM_REF_GENOTYPE_INDEX = 0;

    private final double refPseudocount;
    private final double snpPseudocount;
    private final double indelPseudocount;
    private final int defaultPloidy;


    public AlleleFrequencyCalculator(final double refPseudocount, final double snpPseudocount, final double indelPseudocount, final int defaultPloidy) {
        this.refPseudocount = refPseudocount;
        this.snpPseudocount = snpPseudocount;
        this.indelPseudocount = indelPseudocount;
        this.defaultPloidy = defaultPloidy;
    }

    public AFCalculationResult getLog10PNonRef(final VariantContext vc) {
        // maxAltAlleles is not used by getLog10PNonRef, so don't worry about the 0
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

        double[] alleleCounts = new double[numAlleles];
        final double flatLog10AlleleFrequency = -MathUtils.log10(numAlleles); // log10(1/numAlleles)
        double[] log10AlleleFrequencies = new IndexRange(0, numAlleles).mapToDouble(n -> flatLog10AlleleFrequency);

        for (double alleleCountsMaximumDifference = Double.POSITIVE_INFINITY; alleleCountsMaximumDifference > THRESHOLD_FOR_ALLELE_COUNT_CONVERGENCE; ) {
            final double[] newAlleleCounts = effectiveAlleleCounts(vc, log10AlleleFrequencies);
            alleleCountsMaximumDifference = Arrays.stream(MathArrays.ebeSubtract(alleleCounts, newAlleleCounts)).map(Math::abs).max().getAsDouble();
            alleleCounts = newAlleleCounts;
            final double[] posteriorPseudocounts = MathArrays.ebeAdd(priorPseudocounts, alleleCounts);

            // first iteration uses flat prior in order to avoid local minimum where the prior + no pseudocounts gives such a low
            // effective allele frequency that it overwhelms the genotype likelihood of a real variant
            // basically, we want a chance to get non-zero pseudocounts before using a prior that's biased against a variant
            log10AlleleFrequencies = new Dirichlet(posteriorPseudocounts).log10MeanWeights();
        }

        double[] log10POfZeroCountsByAllele = new double[numAlleles];
        double log10PNoVariant = 0;

        final boolean spanningDeletionPresent = alleles.contains(Allele.SPAN_DEL);
        final Map<Integer, int[]> nonVariantIndicesByPloidy = new Int2ObjectArrayMap<>();
        for (final Genotype g : vc.getGenotypes()) {
            if (!g.hasLikelihoods()) {
                continue;
            }
            final int ploidy = g.getPloidy() == 0 ? defaultPloidy : g.getPloidy();
            final GenotypeLikelihoodCalculator glCalc = GL_CALCS.getInstance(ploidy, numAlleles);

            final double[] log10GenotypePosteriors = log10NormalizedGenotypePosteriors(g, glCalc, log10AlleleFrequencies);

            //the total probability
            if (!spanningDeletionPresent) {
                log10PNoVariant += log10GenotypePosteriors[HOM_REF_GENOTYPE_INDEX];
            } else {
                nonVariantIndicesByPloidy.computeIfAbsent(ploidy, p -> genotypeIndicesWithOnlyRefAndSpanDel(p, alleles));
                final int[] nonVariantIndices = nonVariantIndicesByPloidy.get(ploidy);
                final double[] nonVariantLog10Posteriors = MathUtils.applyToArray(nonVariantIndices, n -> log10GenotypePosteriors[n]);
                // when the only alt allele is the spanning deletion the probability that the site is non-variant
                // may be so close to 1 that finite precision error in log10SumLog10 yields a positive value,
                // which is bogus.  Thus we cap it at 0.
                log10PNoVariant += Math.min(0,MathUtils.log10SumLog10(nonVariantLog10Posteriors));
            }

            // per allele non-log space probabilities of zero counts for this sample
            // for each allele calculate the total probability of genotypes containing at least one copy of the allele
            final double[] log10ProbabilityOfNonZeroAltAlleles = new double[numAlleles];
            Arrays.fill(log10ProbabilityOfNonZeroAltAlleles, Double.NEGATIVE_INFINITY);

            for (int genotype = 0; genotype < glCalc.genotypeCount(); genotype++) {
                final double log10GenotypePosterior = log10GenotypePosteriors[genotype];
                glCalc.genotypeAlleleCountsAt(genotype).forEachAlleleIndexAndCount((alleleIndex, count) ->
                        log10ProbabilityOfNonZeroAltAlleles[alleleIndex] =
                                MathUtils.log10SumLog10(log10ProbabilityOfNonZeroAltAlleles[alleleIndex], log10GenotypePosterior));
            }

            for (int allele = 0; allele < numAlleles; allele++) {
                // if prob of non hom ref == 1 up to numerical precision, short-circuit to avoid NaN
                if (log10ProbabilityOfNonZeroAltAlleles[allele] >= 0) {
                    log10POfZeroCountsByAllele[allele] = Double.NEGATIVE_INFINITY;
                } else {
                    log10POfZeroCountsByAllele[allele] += MathUtils.log10OneMinusPow10(log10ProbabilityOfNonZeroAltAlleles[allele]);
                }
            }
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
        final double[] log10PosteriorOfNoVariantYesVariant = {log10PNoVariant, MathUtils.log10OneMinusPow10(log10PNoVariant)};

        return new AFCalculationResult(integerAltAlleleCounts, alleles, log10PosteriorOfNoVariantYesVariant, dummyFlatPrior, log10PRefByAllele);
    }

    // effectiveAlleleCounts[allele a] = SUM_{genotypes g} (posterior_probability(g) * num_copies of a in g), which we denote as SUM [n_g p_g]
    // for numerical stability we will do this in log space:
    // count = SUM 10^(log (n_g p_g)) = SUM 10^(log n_g + log p_g)
    // thanks to the log-sum-exp trick this lets us work with log posteriors alone
    private double[] effectiveAlleleCounts(final VariantContext vc, final double[] log10AlleleFrequencies) {
        final int numAlleles = vc.getNAlleles();
        Utils.validateArg(numAlleles == log10AlleleFrequencies.length, "number of alleles inconsistent");
        final double[] log10Result = new double[numAlleles];
        Arrays.fill(log10Result, Double.NEGATIVE_INFINITY);
        for (final Genotype g : vc.getGenotypes()) {
            if (!g.hasLikelihoods()) {
                continue;
            }
            final GenotypeLikelihoodCalculator glCalc = GL_CALCS.getInstance(g.getPloidy(), numAlleles);

            final double[] log10GenotypePosteriors = log10NormalizedGenotypePosteriors(g, glCalc, log10AlleleFrequencies);

            new IndexRange(0, glCalc.genotypeCount()).forEach(genotypeIndex ->
                glCalc.genotypeAlleleCountsAt(genotypeIndex).forEachAlleleIndexAndCount((alleleIndex, count) ->
                        log10Result[alleleIndex] = MathUtils.log10SumLog10(log10Result[alleleIndex], log10GenotypePosteriors[genotypeIndex] + MathUtils.log10(count))));
        }
        return MathUtils.applyToArrayInPlace(log10Result, x -> Math.pow(10.0, x));
    }

    private static double[] log10NormalizedGenotypePosteriors(final Genotype g, final GenotypeLikelihoodCalculator glCalc, final double[] log10AlleleFrequencies) {
        final double[] log10Likelihoods = g.getLikelihoods().getAsVector();
        final double[] log10Posteriors = new IndexRange(0, glCalc.genotypeCount()).mapToDouble(genotypeIndex -> {
            final GenotypeAlleleCounts gac = glCalc.genotypeAlleleCountsAt(genotypeIndex);
            return gac.log10CombinationCount() + log10Likelihoods[genotypeIndex]
                    + gac.sumOverAlleleIndicesAndCounts((index, count) -> count * log10AlleleFrequencies[index]);
        });
        return MathUtils.normalizeLog10(log10Posteriors);
    }

    private static int[] genotypeIndicesWithOnlyRefAndSpanDel(final int ploidy, final List<Allele> alleles) {
        final GenotypeLikelihoodCalculator glCalc = GL_CALCS.getInstance(ploidy, alleles.size());
        final boolean spanningDeletionPresent = alleles.contains(Allele.SPAN_DEL);
        if (!spanningDeletionPresent) {
            return new int[] {HOM_REF_GENOTYPE_INDEX};
        } else {
            final int spanDelIndex = alleles.indexOf(Allele.SPAN_DEL);
            // allele counts are in the GenotypeLikelihoodCalculator format of {ref index, ref count, span del index, span del count}
            return new IndexRange(0, ploidy).mapToInteger(n -> glCalc.alleleCountsToIndex(new int[]{0, ploidy - n, spanDelIndex, n}));
        }
    }

    @Override   //Note: unused
    protected AFCalculationResult getResultFromFinalState(final VariantContext vc, final double[] priors, final StateTracker st) { return null; }

    @Override//Note: unused
    protected AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy,
                                                               final double[] priors, final StateTracker st) { return null; }

    @Override   //Note: unused
    protected StateTracker getStateTracker(final boolean reset, final int maximumAlternativeAlleleCount) { return null; }

}
