package org.broadinstitute.hellbender.tools.funcotator.filtrationRules;


import htsjdk.variant.variantcontext.VariantContext;

public abstract class TwoPassFuncotationFilter extends FuncotationFilter {
    TwoPassFuncotationFilter(String filterName) {
        super(filterName);
    }

    public abstract void firstPassApply(final VariantContext variant);

    public abstract void afterFirstPass();
}
