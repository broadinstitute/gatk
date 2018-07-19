package org.broadinstitute.hellbender.tools.funcotator.filtrationRules;

import java.util.Collections;
import java.util.List;

/**
 * {@link FuncotationFilter} matching variants which:
 * <ul>
 *     <li>Have been flagged by LMM as important for loss of function.</li>
 * </ul>
 */
public class LmmFilter extends FuncotationFilter {

    /**
     * Value to include in the {@value org.broadinstitute.hellbender.tools.funcotator.FilterFuncotationsConstants#CLINSIG_INFO_KEY}
     * INFO annotation of variants matching this rule.
     */
    public static final String CLINSIG_INFO_VALUE = "LMM";

    /**
     * Funcotation which will contain "true" for variants which LMM has marked as important.
     */
    private static final String LMM_FLAGGED = "LMMKnown_LMM_FLAGGED";

    public LmmFilter() {
        super(CLINSIG_INFO_VALUE);
    }

    @Override
    List<FuncotationFiltrationRule> getRules() {
        return Collections.singletonList(funcotations -> Boolean.valueOf(funcotations.getOrDefault(LMM_FLAGGED, "false")));
    }
}
