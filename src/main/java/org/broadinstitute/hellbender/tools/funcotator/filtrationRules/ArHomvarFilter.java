package org.broadinstitute.hellbender.tools.funcotator.filtrationRules;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class ArHomvarFilter extends FuncotationFilter {
    public ArHomvarFilter() {
        super(AutosomalRecessiveConstants.AR_INFO_VALUE);
    }

    @Override
    List<FuncotationFiltrationRule> getRules() {
        return Collections.singletonList(this::arHomvarRule);
    }

    private boolean arHomvarRule(Set<Map.Entry<String, String>> funcotations, VariantContext variantContext) {
        // Is this gene part of the list of genes we care about?
        boolean isInterestingGene = funcotations.stream()
                .anyMatch(funcotation ->
                        funcotation.getKey().equals("Gencode_27_hugoSymbol")
                        && AutosomalRecessiveConstants.AUTOSOMAL_RECESSIVE_GENES.contains(funcotation.getValue())
                );
        // If so, is this a homvar?
        if (isInterestingGene) {
            return variantContext.getHomVarCount() > 0;
        }

        return false;
    }
}
