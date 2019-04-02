package org.broadinstitute.hellbender.tools.funcotator.filtrationRules;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.funcotator.FilterFuncotations;
import org.broadinstitute.hellbender.tools.funcotator.FilterFuncotationsUtils;
import org.broadinstitute.hellbender.tools.funcotator.FuncotationMap;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

public class ArHetvarFilter extends TwoPassFuncotationFilter {
    /**
     * Funcotation which will contain the gene on which this variant lives, determined by Funcotator.
     *
     * Varies based on gencode version.
     */
    private final String gene;

    private final String annotationTranscript;
    private final List<VariantContext> arCompoundHetVariants = new ArrayList<>();
    private final Map<String, List<VariantContext>> arHetVariantsByGene = new HashMap<>();
    private final String[] funcotationKeys;
    private boolean firstPassApplied;
    private boolean afterFirstPassApplied;


    @Override
    List<FuncotationFiltrationRule> getRules() {
        return Collections.singletonList(this::arHetvarRule);
    }

    public ArHetvarFilter(FilterFuncotations.Reference reference, final String[] funcotationKeys) {
        super(AutosomalRecessiveConstants.AR_INFO_VALUE);
        this.gene = "Gencode_" + reference.getGencodeVersion() + "_hugoSymbol";
        this.annotationTranscript = "Gencode_" + reference.getGencodeVersion() + "_annotationTranscript";
        this.funcotationKeys = funcotationKeys;
    }

    @Override
    public void firstPassApply(final VariantContext variant) {
        firstPassApplied = true;
        buildArHetByGene(variant);
    }

    @Override
    public void afterFirstPass() {
        if (!firstPassApplied) {
            throw new GATKException("firstPassApply should be called before afterFirstPass");
        }
        afterFirstPassApplied = true;
        arHetVariantsByGene.keySet().forEach(gene -> {
            if(arHetVariantsByGene.get(gene).size() > 1) {
                arCompoundHetVariants.addAll(arHetVariantsByGene.get(gene));
            }
        });
    }

    private boolean arHetvarRule(Set<Map.Entry<String, String>> funcotations, VariantContext variant) {
        if (!firstPassApplied) {
            throw new GATKException("firstPassApply should be called before this rule is applied");
        }
        else if (!afterFirstPassApplied) {
            throw new GATKException("afterFirstPassApplied should be called before this rule is applied");
        }
        return arCompoundHetVariants.stream().anyMatch(hetVariant -> variantContextsMatch(hetVariant, variant));
    }

    private void buildArHetByGene(final VariantContext variant) {
        final Map<Allele, FuncotationMap> funcs = FuncotatorUtils.createAlleleToFuncotationMapFromFuncotationVcfAttribute(
                funcotationKeys, variant, annotationTranscript, "FAKE_SOURCE");

        funcs.values().forEach(funcotationMap -> {
            FilterFuncotationsUtils.getTranscriptFuncotations(funcotationMap).forEach(funcotations -> {
                Optional<Map.Entry<String, String>> maybeGeneFuncotation = funcotations.stream().filter(funcotation -> funcotation.getKey().equals(gene)).findFirst();
                if (maybeGeneFuncotation.isPresent()) {
                    String gene = maybeGeneFuncotation.get().getValue();
                    if (AutosomalRecessiveConstants.AUTOSOMAL_RECESSIVE_GENES.contains(gene) && variant.getHetCount() > 0) {
                        if(arHetVariantsByGene.containsKey(gene)) {
                            arHetVariantsByGene.get(gene).add(variant);
                        }
                        else {
                            ArrayList<VariantContext> variants = new ArrayList<>();
                            variants.add(variant);
                            arHetVariantsByGene.put(gene, variants);
                        }
                    }
                }
            });
        });
    }

    // We know these VariantContexts come from the same list of variants, so we should only need to check these things
    // instead of these things plus attributes, filters, and qual scores.
    private boolean variantContextsMatch(VariantContext v1, VariantContext v2) {
        return v1.getContig().equals(v2.getContig())
                && v1.getStart() == v2.getStart()
                && v1.getEnd() == v2.getEnd()
                && v1.getReference() == v2.getReference()
                && v1.getAlternateAlleles().size() == v2.getAlternateAlleles().size()
                && v1.getAlternateAlleles().containsAll(v2.getAlternateAlleles());
    }
}
