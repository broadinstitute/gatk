package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.Serializable;

/**
 * VariantContextVariantAdapter wraps the existing htsjdk VariantContext class so it can be
 * used with the {@link GATKVariant} API.
 */
public class VariantContextVariantAdapter implements GATKVariant, Serializable {
    private static final long serialVersionUID = 1L;

    private final VariantContext variantContext;

    public VariantContextVariantAdapter(VariantContext vc) {
        this.variantContext = vc;
    }

    public static GATKVariant sparkVariantAdapter(VariantContext vc) {
        return new MinimalVariant(new SimpleInterval(vc.getContig(),vc.getStart(),vc.getEnd()), vc.isSNP(), vc.isIndel());
    }

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
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        VariantContextVariantAdapter that = (VariantContextVariantAdapter) o;

        // VariantContext doesn't define equality, so we have to.
        if (!getContig().equals(that.getContig())) {
            return false;
        }
        if (getStart() != that.getStart()) {
            return false;
        }
        if (getEnd() != that.getEnd()) {
            return false;
        }
        if (isSnp() != that.isSnp()) {
            return false;
        }
        return isIndel() != that.isIndel();
    }

    @Override
    public int hashCode() {
        return variantContext.hashCode();
    }

    @Override
    public String toString() {
        return String.format("VariantContextVariantAdapter -- interval(%s:%d-%d), snp(%b), indel(%b)",
                getContig(), getStart(), getEnd(), isSnp(), isIndel());
    }
}
