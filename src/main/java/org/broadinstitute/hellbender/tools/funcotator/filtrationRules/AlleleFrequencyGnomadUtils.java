package org.broadinstitute.hellbender.tools.funcotator.filtrationRules;


import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

/**
 * Allele frequency calculations for the Gnomad dataset
 */
class AlleleFrequencyGnomadUtils {
    /**
     * Sub-population suffixes used within gnomAD. Used for calculating max MAF.
     */
    private enum GnomadSubpopSuffixes {
        afr, // African / African-American
        afr_female,
        afr_male,
        amr, // Latino
        amr_female,
        amr_male,
        asj, // Ashkenazi Jewish
        asj_female,
        asj_male,
        eas, // East Asian
        eas_female,
        eas_jpn,
        eas_kor,
        eas_male,
        eas_oea,
        female,
        fin, // Finnish
        fin_female,
        fin_male,
        male,
        nfe, // Non-Finnish European
        nfe_bgr,
        nfe_est,
        nfe_female,
        nfe_male,
        nfe_nwe,
        nfe_onf,
        nfe_seu,
        nfe_swe,
        oth, // Other
        oth_female,
        oth_male,
        popmax,
        raw,
        sas, // South Asian
        sas_female,
        sas_male
    }

    /**
     * Prefix for gnomAD genome allele-frequency Funcotations.
     */
    private enum GnomadDataset {
        gnomAD_genome, gnomAD_exome;

        // Prefix for gnomAD allele-frequency Funcotations.
        private String alleleFrequencyPrefix = name() + "_AF_";

        // Name of annotations for whether gnomAD records for a particular data type were filtered out.
        private String filterAnnotation = name() + "_FILTER";
    }

    /**
     * Calculate the max MAF across all gnomAD sub-populations from the given Funcotations.
     * If a sub-population has an allele number of zero, it will be assigned a MAF of zero.
     */
    static double getMaxMinorAlleleFreq(final Map<String, List<String>> funcotations) {
        return getSubpopAlleleFrequencies(funcotations)
                .filter(funcotations::containsKey)
                .flatMapToDouble(subpop -> funcotations.get(subpop).stream().mapToDouble(Double::valueOf))
                .max()
                .orElse(0d);
    }

    static boolean allFrequenciesFiltered(final Map<String, List<String>> funcotations) {
        return datasetsPresent(funcotations).count() == 0;
    }

    private static Stream<String> getSubpopAlleleFrequencies(final Map<String, List<String>> funcotations) {
        return datasetsPresent(funcotations)
                .flatMap(dataset -> Arrays.stream(GnomadSubpopSuffixes.values()).map(suffix -> dataset.alleleFrequencyPrefix + suffix.name()));
    }

    private static Stream<GnomadDataset> datasetsPresent(final Map<String, List<String>> funcotations) {
        return Arrays.stream(GnomadDataset.values())
                .filter(dataset ->
                        !funcotations.containsKey(dataset.filterAnnotation) ||
                                (funcotations.containsKey(dataset.filterAnnotation) && funcotations.get(dataset.filterAnnotation).stream().anyMatch(a -> a.equals("PASS"))));
    }
}
