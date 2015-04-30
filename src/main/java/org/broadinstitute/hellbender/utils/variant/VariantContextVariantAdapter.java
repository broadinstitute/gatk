package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.variantcontext.VariantContext;

import java.io.Serializable;

/**
 * VariantContextVariantAdapter wraps the existing htsjdk VariantContext class so it can be
 * used with the Variant API.
 */
public class VariantContextVariantAdapter implements Variant, Serializable {
    final private VariantContext variantContext;

    public VariantContextVariantAdapter( VariantContext vc ) { this.variantContext = vc; }
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
    @Override
    public boolean equals(Object o) {
        if (!(o instanceof Variant)) {
            return false;
        }

        return ((Variant) o).getContig().equals(getContig()) &&
                ((Variant) o).getStart() == getStart() &&
                ((Variant) o).getEnd() == getEnd() &&
                ((Variant) o).isSnp() == isSnp() &&
                ((Variant) o).isIndel() == isIndel();
    }
    @Override
    public int hashCode() {
        // Per Effective Java Item 9, use the pattern c' = c + 31*x, where x is the next feature.
        // This rule is applied for every field.
        return getContig().hashCode() + 31 *
                (getStart() + 31 * (getEnd() + 31 * (Boolean.hashCode(isSnp())) + 31 * Boolean.hashCode(isIndel())));
    }
}