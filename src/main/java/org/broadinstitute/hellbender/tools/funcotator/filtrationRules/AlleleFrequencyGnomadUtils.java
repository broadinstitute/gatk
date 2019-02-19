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
        _afr, // African / African-American
        _afr_female,
        _afr_male,
        _amr, // Latino
        _amr_female,
        _amr_male,
        _asj, // Ashkenazi Jewish
        _asj_female,
        _asj_male,
        _eas, // East Asian
        _eas_female,
        _eas_jpn,
        _eas_kor,
        _eas_male,
        _eas_oea,
        _female,
        _fin, // Finnish
        _fin_female,
        _fin_male,
        _male,
        _nfe, // Non-Finnish European
        _nfe_bgr,
        _nfe_est,
        _nfe_female,
        _nfe_male,
        _nfe_nwe,
        _nfe_onf,
        _nfe_seu,
        _nfe_swe,
        _oth, // Other
        _oth_female,
        _oth_male,
        _popmax,
        _raw,
        _sas, // South Asian
        _sas_female,
        _sas_male
    }

    /**
     * Prefix for gnomAD genome allele-frequency Funcotations.
     */
    private static String GNOMAD_GENOME_PREFIX = "gnomAD_genome_AF";

    /**
     * Prefix for gnomAD exome allele-frequency Funcotations.
     */
    private static String GNOMAD_EXOME_PREFIX = "gnomAD_exome_AF";

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
    protected static double getMaxMinorAlleleFreq(final Set<Map.Entry<String, String>> funcotations) {
        List<String> prefixes = Arrays.asList(GNOMAD_EXOME_PREFIX, GNOMAD_GENOME_PREFIX);
        Map<String, String> condensedFuncotations = funcotations.stream().collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
        if (condensedFuncotations.containsKey(GNOMAD_GENOME_FILTER) && !condensedFuncotations.get(GNOMAD_GENOME_FILTER).equals("PASS")) {
            prefixes.remove(GNOMAD_GENOME_PREFIX);
        }
        if (condensedFuncotations.containsKey(GNOMAD_EXOME_FILTER) && !condensedFuncotations.get(GNOMAD_EXOME_FILTER).equals("PASS")) {
            prefixes.remove(GNOMAD_EXOME_PREFIX);
        }

        return Arrays.stream(GnomadSubpopSuffixes.values())
                .flatMap(suffix -> Arrays.stream(prefixes.toArray()).map(prefix -> prefix + suffix.name()))
                .filter(condensedFuncotations::containsKey)
                .map(subpop -> Double.valueOf(condensedFuncotations.get(subpop)))
                .max(Double::compareTo)
                .orElse(0d);
    }
}
