package org.broadinstitute.hellbender.tools.funcotator.filtrationRules;

import org.broadinstitute.hellbender.tools.funcotator.FilterFuncotations;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

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
     * Variant classifications which should be matched by this filter.
     */
    private static final Set<String> CONSTANT_LOF_CLASSIFICATIONS = Stream.of(
            GencodeFuncotation.VariantClassification.FRAME_SHIFT_DEL,
            GencodeFuncotation.VariantClassification.FRAME_SHIFT_INS,
            GencodeFuncotation.VariantClassification.NONSENSE,
            GencodeFuncotation.VariantClassification.START_CODON_DEL,
            GencodeFuncotation.VariantClassification.SPLICE_SITE
    ).map(GencodeFuncotation.VariantClassification::toString).collect(Collectors.toSet());

    /**
     * Maximum MAF a variant can have in ExAC to pass this rule.
     */
    private static final double LOF_MAX_MAF = 0.01;

    /**
     * Funcotation which will contain the variant classification determined by Funcotator.
     *
     * Varies based on gencode version.
     */
    private final String classificationFuncotation;

    public LofFilter(final FilterFuncotations.Reference ref) {
        super(CLINSIG_INFO_VALUE);
        this.classificationFuncotation = "Gencode_" + ref.getGencodeVersion() + "_variantClassification";
    }

    @Override
    List<FuncotationFiltrationRule> getRules() {
        return Arrays.asList(
                funcotations -> CONSTANT_LOF_CLASSIFICATIONS.contains(funcotations.getOrDefault(classificationFuncotation, "")),
                funcotations -> funcotations.getOrDefault(LOF_GENE_FUNCOTATION, "").equals("YES"),
                FilterFuncotationsExacUtils.buildExacMaxMafRule(LOF_MAX_MAF));
    }
}
