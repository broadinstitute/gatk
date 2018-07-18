package org.broadinstitute.hellbender.tools.funcotator.filtrationRules;

import java.util.Arrays;
import java.util.Map;

public class FilterFuncotationsExacUtils {
    /**
     * Sub-population suffixes used within ExAC. Used for calculating max MAF.
     */
    private enum ExacSubPopulation {
        AFR, AMR, EAS, FIN, NFE, OTH, SAS
    }

    /**
     * Prefix for allele-count Funcotations for each ExAC sub-population.
     */
    private static String EXAC_ALLELE_COUNT_PREFIX = "ExAC_AC_";

    /**
     * Prefix for allele-number Funcotations for each ExAC sub-population.
     */
    private static String EXAC_ALLELE_NUMBER_PREFIX = "ExAC_AN_";

    /**
     * Build a {@link FuncotationFiltrationRule} matching variants with a MAF less than
     * the given value across all sub-populations of ExAC.
     */
    public static FuncotationFiltrationRule buildExacMaxMafRule(final double maxMaf) {
        return ((variant, prunedTranscriptFuncotations) -> getMaxMinorAlleleFreq(prunedTranscriptFuncotations) <= maxMaf);
    }

    /**
     * Calculate the max MAF across all ExAC sub-populations from the given Funcotations.
     *
     * If a sub-population has an allele number of zero, it will be assigned a MAF of zero.
     */
    private static double getMaxMinorAlleleFreq(final Map<String, String> funcotations) {
        return Arrays.stream(ExacSubPopulation.values())
                .filter(subpop -> funcotations.containsKey(EXAC_ALLELE_COUNT_PREFIX + subpop.name()))
                .map(subpop -> {
                    final Double ac = Double.valueOf(funcotations.get(EXAC_ALLELE_COUNT_PREFIX + subpop.name()));
                    final Integer an = Integer.valueOf(funcotations.get(EXAC_ALLELE_NUMBER_PREFIX + subpop.name()));

                    if (an == 0) {
                        // If a variant has never been seen in ExAC, report it as 0% MAF.
                        return 0d;
                    } else {
                        return ac / an;
                    }
                })
                .max(Double::compareTo)
                .orElse(0d);
    }
}
