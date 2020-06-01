package org.broadinstitute.hellbender.tools.haplotypecaller;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAlleleCounts;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.utils.pairhmm.DragstrParams;

import java.util.Arrays;
import java.util.List;

public class SimpleGenotypePriorCalculator implements GenotypePriorCalculator {



    private enum AlleleType {
        REF, SNP, INDEL, OTHER;
    }

    private final double[] hetValues;
    private final double[] homValues;
    private final double[] diffValues;

    private SimpleGenotypePriorCalculator(final double snpHet, final double snpHom,
                                          final double indelHet, final double indelHom,
                                          final double otherHet, final double otherHom) {
        hetValues = new double[4];
        homValues = new double[4];
        diffValues = new double[4];

        // SNPs: // * 1/3 since there is three possible mutations for SNPs.
        hetValues[1] = snpHet - Math.log10(3);
        homValues[1] = snpHom - Math.log10(3);
        // INDELs:
        hetValues[2] = indelHet;
        homValues[2] = indelHom;
        // Others:
        hetValues[3] = otherHet;
        homValues[3] = otherHom;

        for (int i = 0; i < 4; i++) {
            diffValues[i] = homValues[i] - hetValues[i];
        }
    }

    public static SimpleGenotypePriorCalculator givenHetToHomRatio(final double snpHet, final double indelHet,
                                                                   final double otherHet, final double het_hom_ratio) {
        final double log10Ratio = Math.log10(het_hom_ratio);
        return new SimpleGenotypePriorCalculator(snpHet, snpHet - log10Ratio,
                                                  indelHet, indelHet - log10Ratio,
                                                  otherHet, otherHet - log10Ratio);
    }

    public static SimpleGenotypePriorCalculator givenDragstrParams(final DragstrParams dragstrParams, final int period,
                                                                   final int repeats, final double snpHeterozygosity,
                                                                   final double het_hom_ratio) {
        final double snpHet = snpHeterozygosity;
        final double indelHet = -.1 * dragstrParams.api(period, repeats);
        final double otherHet = Math.max(snpHet, indelHet);
        return givenHetToHomRatio(snpHet, indelHet, otherHet, het_hom_ratio);
    }

    public static SimpleGenotypePriorCalculator assumingHW(final double snpHet, final double indelHet) {
        return assumingHW(snpHet, indelHet, Math.max(snpHet, indelHet));
    }

    public static SimpleGenotypePriorCalculator assumingHW(final double snpHet, final double indelHet,
                                                                   final double otherHet) {
        return new SimpleGenotypePriorCalculator(snpHet, snpHet * 2,
                indelHet, indelHet * 2,
                otherHet, otherHet * 2);
    }

    public static SimpleGenotypePriorCalculator assumingHW(GenotypeCalculationArgumentCollection genotypeArgs) {
        return assumingHW(Math.log10(genotypeArgs.snpHeterozygosity),
                          Math.log10(genotypeArgs.indelHeterozygosity));
    }


    @Override
    public double[] getLog10Priors(final GenotypeLikelihoodCalculator lkCalculator, final List<Allele> alleles) {
        final int[] alleleTypes = calculateAlleleTypes(alleles);
        final int numberOfGenotypes = lkCalculator.genotypeCount();
        final double[] result = new double[numberOfGenotypes];
        // implied = result[0] = 0.0;
        for (int g = 1; g < numberOfGenotypes; g++) {
            final GenotypeAlleleCounts gac = lkCalculator.genotypeAlleleCountsAt(g);
            final int numberOfDistictAlleles = gac.distinctAlleleCount();
            double log10Sum = 0;
            for (int a = 0; a < numberOfDistictAlleles; a++) {
                final int idx = gac.alleleIndexAt(a);
                final int cnt = gac.alleleCountAt(a);
                if (cnt == 1) {
                    log10Sum += hetValues[alleleTypes[idx]];
                } else if (cnt == 2) {
                    log10Sum += homValues[alleleTypes[idx]];
                } else { // for plodies over 2 and allele counts over 2 then we use the het/hom ratio for the rest
                    log10Sum += hetValues[alleleTypes[idx]] + diffValues[alleleTypes[idx]] * (cnt - 1);
                }
            }
            result[g] = log10Sum;
        }
        return result;
    }

    private int[] calculateAlleleTypes(final List<Allele> alleles) {
        if (alleles.isEmpty()) {
            throw new IllegalArgumentException("there must be at least one allele (the reference)");
        } else {
            final int[] result = new int[alleles.size()];
            Arrays.fill(result, AlleleType.OTHER.ordinal());
            final Allele refAllele = alleles.get(0);
            if (!refAllele.isReference()) {
                throw new IllegalArgumentException("the first allele in the list must be the reference");
            }
            final int refAlleleLength = refAllele.length();
            result[0] = AlleleType.REF.ordinal();
            for (int i = 1; i < result.length; i++) {
                final Allele allele = alleles.get(i);
                if (allele.isCalled()) {
                    if (!allele.isSymbolic()) {
                        result[i] = (allele.length() == refAlleleLength ? AlleleType.SNP : AlleleType.INDEL).ordinal();
                    } else if (allele.equals(Allele.SV_SIMPLE_INS) || allele.equals(Allele.SV_SIMPLE_DEL)) {
                        result[i] = AlleleType.INDEL.ordinal();
                    }
                }
            }
            return result;
        }
    }
}
