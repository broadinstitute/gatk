package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAlleleCounts;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
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
public final class AlleleFrequencyCalculator {

    private static final GenotypeLikelihoodCalculators GL_CALCS = new GenotypeLikelihoodCalculators();
    private static final double THRESHOLD_FOR_ALLELE_COUNT_CONVERGENCE = 0.1;
    private static final int HOM_REF_GENOTYPE_INDEX = 0;

    private final double refPseudocount;
    private final double snpPseudocount;
    private final double indelPseudocount;
    private final int defaultPloidy;

    public AlleleFrequencyCalculator(final double refPseudocount, final double snpPseudocount,
                                     final double indelPseudocount, final int defaultPloidy) {
        this.refPseudocount = refPseudocount;
        this.snpPseudocount = snpPseudocount;
        this.indelPseudocount = indelPseudocount;
        this.defaultPloidy = defaultPloidy;
    }

    public static AlleleFrequencyCalculator makeCalculator(GenotypeCalculationArgumentCollection genotypeArgs) {
        final double refPseudocount = genotypeArgs.snpHeterozygosity / Math.pow(genotypeArgs.heterozygosityStandardDeviation,2);
        final double snpPseudocount = genotypeArgs.snpHeterozygosity * refPseudocount;
        final double indelPseudocount = genotypeArgs.indelHeterozygosity * refPseudocount;
        return new AlleleFrequencyCalculator(refPseudocount, snpPseudocount, indelPseudocount, genotypeArgs.samplePloidy);
    }

    public static AlleleFrequencyCalculator makeCalculator(final DragstrParams dragstrParms, final int period,
                                                           final int repeats, final int ploidy,
                                                           final double snpHeterozygosity, final double scale) {
        final double api = dragstrParms.api(period, repeats);
        final double log10IndelFreq = api * -.1;
        final double log10RefFreq = MathUtils.log10OneMinusPow10(log10IndelFreq);
        final double log10SnpFreq = log10RefFreq + Math.log10(snpHeterozygosity);
        final double log10Sum = MathUtils.log10SumLog10(log10RefFreq, log10IndelFreq, log10SnpFreq);
        final double log10ScaleUp = Math.log10(scale) - log10Sum;
        final double refPseudoCount = Math.pow(10, log10ScaleUp + log10RefFreq);
        final double indelPseudoCount = Math.pow(10, log10ScaleUp + log10IndelFreq);
        final double snpPseudoCount = Math.pow(10, log10ScaleUp + log10SnpFreq);
        return new AlleleFrequencyCalculator(refPseudoCount, snpPseudoCount, indelPseudoCount, ploidy);
    }

    /**
     *
     * @param g must have likelihoods or (if approximateHomRefsFromGQ is true) GQ
     * @param glCalc
     * @param log10AlleleFrequencies
     * @return
     */
    private static double[] log10NormalizedGenotypePosteriors(final Genotype g, final GenotypeLikelihoodCalculator glCalc, final double[] log10AlleleFrequencies) {
        final double[] log10Likelihoods;
        if (g.hasLikelihoods()) {
            log10Likelihoods = g.getLikelihoods().getAsVector();
        } else if ( g.isHomRef() || g.isNoCall()) {
            if (g.hasGQ()) {
                int[] perSampleIndexesOfRelevantAlleles = new int[log10AlleleFrequencies.length];
                for (int i = 1; i < perSampleIndexesOfRelevantAlleles.length; i++) {
                    perSampleIndexesOfRelevantAlleles[i] = 1;
                }
                final int gq = g.getGQ();
                final int ploidy = g.getPloidy();
                final int[] approxLikelihoods = {0, gq, 15*gq};
                final int[] genotypeIndexMapByPloidy = GL_CALCS.getInstance(ploidy, log10AlleleFrequencies.length).genotypeIndexMap(perSampleIndexesOfRelevantAlleles, GL_CALCS); //probably horribly slow
                final int[] PLs = new int[genotypeIndexMapByPloidy.length];
                for (int i = 0; i < PLs.length; i++) {
                        PLs[i] = approxLikelihoods[genotypeIndexMapByPloidy[i]];
                    }
                log10Likelihoods = GenotypeLikelihoods.fromPLs(PLs).getAsVector();
            } else {
                throw new IllegalStateException("Genotype " + g + " does not contain GQ necessary to calculate posteriors.");
            }
        } else {
            throw new IllegalStateException("Genotype " + g + " does not contain likelihoods necessary to calculate posteriors.");
        }
        final double[] log10Posteriors = new IndexRange(0, glCalc.genotypeCount()).mapToDouble(genotypeIndex -> {
            final GenotypeAlleleCounts gac = glCalc.genotypeAlleleCountsAt(genotypeIndex);
            return gac.log10CombinationCount() + log10Likelihoods[genotypeIndex]
                    + gac.sumOverAlleleIndicesAndCounts((index, count) -> count * log10AlleleFrequencies[index]);
        });
        return MathUtils.normalizeLog10(log10Posteriors);
    }

    private static double[] approximateHomRefPLsFromGQ(final int genotypeCount, final int gq, final GenotypeLikelihoodCalculator glCalc) {
        final int[] newPLs = new int[genotypeCount];
        for (int i = 0; i < genotypeCount; i++) {
            final GenotypeAlleleCounts gac = glCalc.genotypeAlleleCountsAt(i);
            final int refCount = gac.alleleCountFor(0);
                newPLs[i] = (glCalc.ploidy() - refCount) * gq;
        }
        return GenotypeLikelihoods.fromPLs(newPLs).getAsVector();
    }

    private static int[] genotypeIndicesWithOnlyRefAndSpanDel(final int ploidy, final List<Allele> alleles) {
        final GenotypeLikelihoodCalculator glCalc = GL_CALCS.getInstance(ploidy, alleles.size());
        final boolean spanningDeletionPresent = alleles.contains(Allele.SPAN_DEL);
        if (!spanningDeletionPresent) {
            return new int[] {HOM_REF_GENOTYPE_INDEX};
        } else {
            final int spanDelIndex = alleles.indexOf(Allele.SPAN_DEL);
            // allele counts are in the GenotypeLikelihoodCalculator format of {ref index, ref count, span del index, span del count}
            return new IndexRange(0, ploidy + 1).mapToInteger(n -> glCalc.alleleCountsToIndex(new int[]{0, ploidy - n, spanDelIndex, n}));
        }
    }

    public int getPloidy() {
        return defaultPloidy;
    }

    public AFCalculationResult calculate(final VariantContext vc) {
        // maxAltAlleles is not used by getLog10PNonRef, so don't worry about the 0
        return calculate(vc, defaultPloidy);
    }

    /**
     * Compute the probability of the alleles segregating given the genotype likelihoods of the samples in vc
     *
     * @param vc the VariantContext holding the alleles and sample information.  The VariantContext
     *           must have at least 1 alternative allele
     * @return result (for programming convenience)
     */
    public AFCalculationResult calculate(final VariantContext vc, final int defaultPloidy) {
        Utils.nonNull(vc, "VariantContext cannot be null");
        final int numAlleles = vc.getNAlleles();
        final List<Allele> alleles = vc.getAlleles();
        Utils.validateArg( numAlleles > 1, () -> "VariantContext has only a single reference allele, but getLog10PNonRef requires at least one at all " + vc);

        final double[] priorPseudocounts = alleles.stream()
                .mapToDouble(a -> a.isReference() ? refPseudocount : (a.length() == vc.getReference().length() ? snpPseudocount : indelPseudocount)).toArray();

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

        // re-usable buffers of the log10 genotype posteriors of genotypes missing each allele
        final List<DoubleArrayList> log10AbsentPosteriors = IntStream.range(0,numAlleles).mapToObj(n -> new DoubleArrayList()).collect(Collectors.toList());
        for (final Genotype g : vc.getGenotypes()) {
            if (!g.hasLikelihoods() && !g.hasGQ() && (!g.getAlleles().stream().anyMatch(a -> a.isCalled() && a.isNonReference() && !a.isSymbolic()))) {
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

            // if the VC is biallelic the allele-specific qual equals the variant qual
            if (numAlleles == 2 && !spanningDeletionPresent) {
                continue;
            }

            // for each allele, we collect the log10 probabilities of genotypes in which the allele is absent, then add (in log space)
            // to get the log10 probability that the allele is absent in this sample
            log10AbsentPosteriors.forEach(DoubleArrayList::clear);  // clear the buffers.  Note that this is O(1) due to the primitive backing array
            for (int genotype = 0; genotype < glCalc.genotypeCount(); genotype++) {
                final double log10GenotypePosterior = log10GenotypePosteriors[genotype];
                glCalc.genotypeAlleleCountsAt(genotype).forEachAbsentAlleleIndex(a -> log10AbsentPosteriors.get(a).add(log10GenotypePosterior), numAlleles);
            }

            final double[] log10PNoAllele = log10AbsentPosteriors.stream()
                    .mapToDouble(buffer -> MathUtils.log10SumLog10(buffer.toDoubleArray()))
                    .map(x -> Math.min(0, x)).toArray();    // if prob of non hom ref > 1 due to finite precision, short-circuit to avoid NaN

            // multiply the cumulative probabilities of alleles being absent, which is addition of logs
            MathUtils.addToArrayInPlace(log10POfZeroCountsByAllele, log10PNoAllele);
        }

        // for biallelic the allele-specific qual equals the variant qual, and we short-circuited the calculation above
        if (numAlleles == 2 && !spanningDeletionPresent) {
            log10POfZeroCountsByAllele[1] = log10PNoVariant;
        }

        // unfortunately AFCalculationResult expects integers for the MLE.  We really should emit the EM no-integer values
        // which are valuable (eg in CombineGVCFs) as the sufficient statistics of the Dirichlet posterior on allele frequencies
        final int[] integerAlleleCounts = Arrays.stream(alleleCounts).mapToInt(x -> (int) Math.round(x)).toArray();
        final int[] integerAltAlleleCounts = Arrays.copyOfRange(integerAlleleCounts, 1, numAlleles);

        //skip the ref allele (index 0)
        final Map<Allele, Double> log10PRefByAllele = IntStream.range(1, numAlleles).boxed()
                .collect(Collectors.toMap(alleles::get, a -> log10POfZeroCountsByAllele[a]));

        return new AFCalculationResult(integerAltAlleleCounts, alleles, log10PNoVariant, log10PRefByAllele);
    }

    /**
     * Calculate the posterior probability that a single biallelic genotype is non-ref
     *
     * The nth genotype (n runs from 0 to the sample ploidy, inclusive) contains n copies of the alt allele
     * @param log10GenotypeLikelihoods
     */
    public double calculateSingleSampleBiallelicNonRefPosterior(final double[] log10GenotypeLikelihoods, final boolean returnZeroIfRefIsMax) {
        Utils.nonNull(log10GenotypeLikelihoods);

        if (returnZeroIfRefIsMax && MathUtils.maxElementIndex(log10GenotypeLikelihoods) == 0) {
            return 0;
        }

        final int ploidy = log10GenotypeLikelihoods.length - 1;

        final double[] log10UnnormalizedPosteriors = new IndexRange(0, ploidy + 1)
                .mapToDouble(n -> log10GenotypeLikelihoods[n] + MathUtils.log10BinomialCoefficient(ploidy, n)
                        + MathUtils.logToLog10(Gamma.logGamma(n + snpPseudocount ) + Gamma.logGamma(ploidy - n + refPseudocount)));

        return (returnZeroIfRefIsMax && MathUtils.maxElementIndex(log10UnnormalizedPosteriors) == 0) ? 0.0 :
                1 - MathUtils.normalizeFromLog10ToLinearSpace(log10UnnormalizedPosteriors)[0];
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
            if (!g.hasLikelihoods() && !g.hasGQ() && (!g.getAlleles().stream().anyMatch(a -> a.isCalled() && a.isNonReference() && !a.isSymbolic()))) {
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
}
