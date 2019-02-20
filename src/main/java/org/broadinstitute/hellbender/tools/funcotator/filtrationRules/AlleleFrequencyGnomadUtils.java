package org.broadinstitute.hellbender.tools.funcotator.filtrationRules;

import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class AlleleFrequencyGnomadUtils {
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
    private static String GNOMAD_GENOME_PREFIX = "gnomAD_genome_AF_";

    /**
     * Prefix for gnomAD exome allele-frequency Funcotations.
     */
    private static String GNOMAD_EXOME_PREFIX = "gnomAD_exome_AF_";

    /**
     * Names of annotations for whether gnomAD records for a particular data type were filtered out.
     */
    private static String GNOMAD_GENOME_FILTER = "gnomAD_genome_FILTER";
    private static String GNOMAD_EXOME_FILTER = "gnomAD_exome_FILTER";

    /**
     * Calculate the max MAF across all gnomAD sub-populations from the given Funcotations.
     *
     * If a sub-population has an allele number of zero, it will be assigned a MAF of zero.
     */
    static double getMaxMinorAlleleFreq(final Map<String, String> funcotations) {
        ArrayList<String> prefixes = getPrefixes(funcotations);
        return Arrays.stream(GnomadSubpopSuffixes.values())
                .flatMap(suffix -> Arrays.stream(prefixes.toArray()).map(prefix -> prefix + suffix.name()))
                .filter(funcotations::containsKey)
                .map(subpop -> Double.valueOf(funcotations.get(subpop)))
                .max(Double::compareTo)
                .orElse(0d);
    }

    static boolean allFrequenciesFiltered(final Map<String, String> funcotations) {
        ArrayList<String> prefixes = getPrefixes(funcotations);
        return prefixes.size() == 0;
    }

    private static ArrayList<String> getPrefixes(Map<String, String> funcotations) {
        ArrayList<String> prefixes = new ArrayList<>();
        prefixes.add(GNOMAD_EXOME_PREFIX);
        prefixes.add(GNOMAD_GENOME_PREFIX);

        if (funcotations.containsKey(GNOMAD_EXOME_FILTER) && !funcotations.get(GNOMAD_EXOME_FILTER).equals("PASS")) {
            prefixes.remove(GNOMAD_EXOME_PREFIX);
        }

        if (funcotations.containsKey(GNOMAD_GENOME_FILTER) && !funcotations.get(GNOMAD_GENOME_FILTER).equals("PASS")) {
            prefixes.remove(GNOMAD_GENOME_PREFIX);
        }
        return prefixes;
    }
}
