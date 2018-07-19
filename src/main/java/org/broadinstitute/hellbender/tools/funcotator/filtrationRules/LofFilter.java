package org.broadinstitute.hellbender.tools.funcotator.filtrationRules;

import org.broadinstitute.hellbender.tools.funcotator.FilterFuncotations;

import java.util.Arrays;
import java.util.List;

/**
 * {@link FuncotationFilter} matching variants which:
 * <ul>
 *     <li>Are classified as FRAME_SHIFT_*, NONSENSE, START_CODON_DEL, or SPLICE_SITE</li>
 *     <li>Occur on a gene where loss of function is a disease mechanism</li>
 *     <li>Have a max MAF of 1% across sub-populations of ExAC</li>
 * </ul>
 */
public class LofFilter extends FuncotationFilter {

    /**
     * Value to include in the {@value org.broadinstitute.hellbender.tools.funcotator.FilterFuncotationsConstants#CLINSIG_INFO_KEY}
     * INFO annotation of variants matching this rule.
     */
    public static final String CLINSIG_INFO_VALUE = "LOF";

    /**
     * Funcotation which will contain "YES" for variants which are important for loss of function.
     */
    private static final String LOF_GENE_FUNCOTATION = "ACMGLMMLof_LOF_Mechanism";

    /**
     * Prefix for frame-shift variant classifications which should be matched by this filter.
     */
    private static final String FRAME_SHIFT_PREFIX = "FRAME_SHIFT_";

    /**
     * Variant classifications which should be matched by this filter.
     */
    private static final List<String> CONSTANT_LOF_CLASSIFICATIONS = Arrays.asList(
            "NONSENSE", "START_CODON_DEL", "SPLICE_SITE");

    /**
     * Funcotation which will contain the variant classification determined by Funcotator.
     *
     * Varies based on gencode version.
     */
    private final String classificationFuncotation;

    public LofFilter(final FilterFuncotations.ReferenceVersion ref) {
        super(CLINSIG_INFO_VALUE);
        this.classificationFuncotation = "Gencode_" + ref.getGencodeVersion() + "_variantClassification";
    }

    @Override
    List<FuncotationFiltrationRule> getRules() {
        return Arrays.asList(
                funcotations -> {
                    final String classification = funcotations.getOrDefault(classificationFuncotation, "");
                    return classification.startsWith(FRAME_SHIFT_PREFIX) || CONSTANT_LOF_CLASSIFICATIONS.contains(classification);
                },
                funcotations -> funcotations.getOrDefault(LOF_GENE_FUNCOTATION, "").equals("YES"),
                FilterFuncotationsExacUtils.buildExacMaxMafRule(0.01));
    }
}
