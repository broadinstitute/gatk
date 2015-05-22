package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.variantcontext.VariantContext;

import java.io.Serializable;

/**
 * VariantContextVariantAdapter wraps the existing htsjdk VariantContext class so it can be
 * used with the Variant API.
 */
public class VariantContextVariantAdapter implements Variant, Serializable {
    private static final long serialVersionUID = 1L;

    final private VariantContext variantContext;

    public VariantContextVariantAdapter(VariantContext vc) { this.variantContext = vc; }
    @Override
    public String getContig() { return variantContext.getContig(); }
    @Override
    public int getStart() { return variantContext.getStart(); }
    @Override
    public int getEnd() { return variantContext.getEnd(); }
    @Override
    public boolean isSnp() { return variantContext.isSNP(); }
    @Override
    public boolean isIndel() { return variantContext.isIndel(); }
    @Override
    public String toString() {
        return String.format("VariantContextVariantAdapter -- interval(%s:%d-%d), snp(%b), indel(%b)",
                getContig(), getStart(), getEnd(), isSnp(), isIndel());
    }
}
