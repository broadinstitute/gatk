package org.broadinstitute.hellbender.tools.funcotator.filtrationRules;

import java.util.Arrays;
import java.util.List;

/**
 * {@link FuncotationFilter} matching variants which:
 * <ul>
 *     <li>Occur on a gene in the American College of Medical Genomics (ACMG)'s list of clinically-significant variants</li>
 *     <li>Have been labeled by ClinVar as pathogenic or likely pathogenic</li>
 *     <li>Have a max MAF of 5% across sub-populations of ExAC</li>
 * </ul>
 */
public class ClinVarFilter extends FuncotationFilter {

    /**
     * Value to include in the {@value org.broadinstitute.hellbender.tools.funcotator.FilterFuncotationsConstants#CLINSIG_INFO_KEY}
     * INFO annotation of variants matching this rule.
     */
    public static final String CLINSIG_INFO_VALUE = "CLINVAR";

    /**
     * Funcotation which will be non-empty for variants which occur on a gene in the ACMG's list.
     *
     * @see <a href="https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/">The gene list</a>
     */
    private static final String ACMG_DISEASE_FUNCOTATION = "ACMG_recommendation_Disease_Name";

    /**
     * Funcotation which contains ClinVar's assessment of a variant's clinical significance.
     *
     * @see <a href="https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig">Valid values for significance</a>
     */
    private static final String CLINVAR_SIGNIFICANCE_FUNCOTATION = "ClinVar_VCF_CLNSIG";

    /**
     * Clinically-significant values to check for within the {@value CLINVAR_SIGNIFICANCE_FUNCOTATION} Funcotation.
     */
    private static final List<String> CLINVAR_SIGNIFICANCE_MATCHING_VALUES = Arrays.asList("Pathogenic", "Likely_pathogenic");

    /**
     * Maximum MAF a variant can have in ExAC to pass this rule.
     */
    private static final double CLINVAR_MAX_MAF = 0.05;

    public ClinVarFilter() {
        super(CLINSIG_INFO_VALUE);
    }

    @Override
    List<FuncotationFiltrationRule> getRules() {
        return Arrays.asList(
                funcotations -> funcotations.containsKey(ACMG_DISEASE_FUNCOTATION),
                funcotations -> {
                    final String significance = funcotations.getOrDefault(CLINVAR_SIGNIFICANCE_FUNCOTATION, "");
                    return CLINVAR_SIGNIFICANCE_MATCHING_VALUES.stream().anyMatch(significance::contains);
                },
                FilterFuncotationsExacUtils.buildExacMaxMafRule(CLINVAR_MAX_MAF));
    }
}
