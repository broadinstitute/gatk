package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAlleleCounts;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.List;

/**
 * Class to compose genotype prior probability calculators.
 * 
 * <p>
 *     Contains a collection of static method to create calculators based on different
 *     assumptions and source of knowledge a priori (e.g. {@link #assumingHW(double, double) assumingHW}
 *     or {@link #givenDragstrParams(DragstrParams, int, int, double, double) givenDragstrParams}).
 * </p>
 *
 * <p>
 *     Such priors are obtained by invoking {@link #getLog10Priors(GenotypeLikelihoodCalculator, List).
 *     This method takes on the list of alleles for that variant, an a reference to the genotype likelihood calculator witch determines the ploidy.
 * </p>
 * assumptions
 */
public final class GenotypePriorCalculator {

    private enum AlleleType {
        REF, SNP, INDEL, OTHER
    }

    private static final int NUMBER_OF_ALLELE_TYPES = AlleleType.values().length;

    // A snp can go to 3 different bases (standard-nucs - 1), so we normalize SNP lks accordingly. Here is the
    // log10 constant used for that:
    private static final double LOG10_SNP_NORMALIZATION_CONSTANT =
            MathUtils.log10(Nucleotide.STANDARD_BASES.size() - 1);

    private final double[] hetValues;
    private final double[] homValues;
    private final double[] diffValues;

    private GenotypePriorCalculator(final double snpHet, final double snpHom,
                                          final double indelHet, final double indelHom,
                                          final double otherHet, final double otherHom) {
        hetValues = new double[NUMBER_OF_ALLELE_TYPES];
        homValues = new double[NUMBER_OF_ALLELE_TYPES];

        // REF
        // by convention ref log10 likelihoods are set to 0.
        // so they are already set.

        // SNPs: normalized for all possible mutations (number of nucs (4) - 1)
        hetValues[AlleleType.SNP.ordinal()] = snpHet - LOG10_SNP_NORMALIZATION_CONSTANT;
        homValues[AlleleType.SNP.ordinal()] = snpHom - LOG10_SNP_NORMALIZATION_CONSTANT;
        // INDELs:
        hetValues[AlleleType.INDEL.ordinal()] = indelHet;
        homValues[AlleleType.INDEL.ordinal()] = indelHom;
        // Others:
        hetValues[AlleleType.OTHER.ordinal()] = otherHet;
        homValues[AlleleType.OTHER.ordinal()] = otherHom;

        diffValues = MathArrays.ebeSubtract(homValues, hetValues);
    }

    /**
     * Calculate priors based on fix heterozygosities (per event type) and het to hom-var prior ratio.
     *
     * @param log10SnpHet snp heterozygosity in log10 scale.
     * @param log10IndelHet indel heterozygosity in log10 scale.
     * @param log10OtherHet heterozygosity for other type of variants in log10 scale.
     * @param hetHomRatio ratio between the het-var and hom-var genotype priors for the same allele in linear scale.
     * @return never {@code null}.
     */
    public static GenotypePriorCalculator givenHetToHomRatio(final double log10SnpHet, final double log10IndelHet,
                                                             final double log10OtherHet, final double hetHomRatio) {
        final double log10Ratio = Math.log10(hetHomRatio);
        return new GenotypePriorCalculator(log10SnpHet, log10SnpHet - log10Ratio,
                log10IndelHet, log10IndelHet - log10Ratio,
                log10OtherHet, log10OtherHet - log10Ratio);
    }

    /**
     * Composes a calculator based on Dragstr Parameters.
     * @param dragstrParams the input DRAGstr parameter values.
     * @param period the period of the repeat of interest.
     * @param repeats the repeat length (number of full repeat units) that composes the STR.
     * @param snpHeterozygosity the snp-heterozygosity.
     * @param hetHomRatio the expected prior ratio between het and homVar  (het over homVar).
     * @return never {@code null}.
     */
    public static GenotypePriorCalculator givenDragstrParams(final DragstrParams dragstrParams, final int period,
                                                             final int repeats, final double snpHeterozygosity,
                                                             final double hetHomRatio) {
        final double indelHet = -.1 * dragstrParams.api(period, repeats);
        final double otherHet = Math.max(snpHeterozygosity, indelHet);
        return givenHetToHomRatio(snpHeterozygosity, indelHet, otherHet, hetHomRatio);
    }

    /**
     * Composes a calculator based on Hardy-Weinberg equilibrium so that only the het-priors
     * are need to calculate the rest.
     * @param snpHet the prior for an SNP alternative allele in log10 scale.
     * @param indelHet the prior for an INDEL alternative allele in log10 scale.
     * @return never {@code null}.
     */
    public static GenotypePriorCalculator assumingHW(final double snpHet, final double indelHet) {
        ParamUtils.isNegativeOrZero(snpHet, "snp-het in log10 scale must be 0 or a negative");
        ParamUtils.isNegativeOrZero(snpHet, "indel-het in log10 scale must be 0 or a negative");
        return assumingHW(snpHet, indelHet, Math.max(snpHet, indelHet));
    }

    private static GenotypePriorCalculator assumingHW(final double snpHet, final double indelHet,
                                                     final double otherHet) {
        return new GenotypePriorCalculator(snpHet, snpHet * 2,
                indelHet, indelHet * 2,
                otherHet, otherHet * 2);
    }

    /**
     * Composes a calculator based on Hardy-Weinberg equilibrium so that only the het-priors
     * are need to calculate the rest.
     * @param genotypeArgs source for the snp and iondel heterozygosity to be used.
     * @return never {@code null}.
     */
    public static GenotypePriorCalculator assumingHW(final GenotypeCalculationArgumentCollection genotypeArgs) {
        return assumingHW(Math.log10(genotypeArgs.snpHeterozygosity),
                Math.log10(genotypeArgs.indelHeterozygosity));
    }

    /**
     * Calculates the priors given the alleles to genetype and a likelihood calculator that determines the ploidy
     * of the sample at that site.
     * @param lkCalculator the input calculator
     * @param alleles the input alleles.
     * @throws IllegalArgumentException if either input is {@code null} or the calculator maximum number of supported alleles is less that the input allele size.
     * @return never {@code null}, the array will have as many positions as necessary to hold the priors of all possible
     * unphased genotypes as per the number of input alleles and the input calculator's ploidy.
     */
    public double[] getLog10Priors(final GenotypeLikelihoodCalculator lkCalculator, final List<Allele> alleles) {
        Utils.nonNull(lkCalculator);
        Utils.nonNull(alleles);
        if (lkCalculator.alleleCount() < alleles.size()) {
            throw new IllegalArgumentException("the number of alleles in the input calculator must be at least as large as the number of alleles in the input list");
        }
        final int[] alleleTypes = calculateAlleleTypes(alleles);
        final int numberOfGenotypes = lkCalculator.genotypeCount();
        final double[] result = new double[numberOfGenotypes];
        // implied = result[0] = 0.0;
        for (int g = 1; g < numberOfGenotypes; g++) {
            final GenotypeAlleleCounts gac = lkCalculator.genotypeAlleleCountsAt(g);
            result[g] = gac.sumOverAlleleIndicesAndCounts((idx, cnt) -> cnt == 2
                    ? homValues[alleleTypes[idx]]
                    : hetValues[alleleTypes[idx]] + diffValues[alleleTypes[idx]] * (cnt - 1));
        }
        return result;
    }

    private int[] calculateAlleleTypes(final List<Allele> alleles) {
        Utils.validateArg(!alleles.isEmpty(), "the input allele list cannot empty");
        final Allele refAllele = alleles.get(0);
        Utils.nonNull(refAllele, "the reference allele cannot be null");
        Utils.validateArg(refAllele.isReference(), "the first allele must be a valid reference");
        final int referenceLength = refAllele.length();
        return alleles.stream()
                .map(allele -> {
                    if (allele.isReference()) {
                        return AlleleType.REF;
                    } else if (allele.isCalled() && !allele.isSymbolic()) {
                        return allele.length() == referenceLength ? AlleleType.SNP : AlleleType.INDEL;
                    } else if (allele.equals(Allele.SV_SIMPLE_INS) || allele.equals(Allele.SV_SIMPLE_DEL)) {
                        throw new IllegalArgumentException("cannot handle symbolic indels: " + allele);
                    } else {
                        return AlleleType.OTHER;
                    }
                })
                .mapToInt(AlleleType::ordinal)
                .toArray();
    }
}
