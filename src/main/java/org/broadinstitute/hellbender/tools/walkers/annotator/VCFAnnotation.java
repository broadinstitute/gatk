package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;

import java.util.function.Function;

public interface VCFAnnotation<T> {

    VCFAnnotationMeta meta();

    default T of(final VariantContext variant, final T defaultValue) {
        if (meta().target() != VCFAnnotationTarget.VARIANT) {
            throw new UnsupportedOperationException("wrong annotation target type");
        } else {
            final Object value = variant.getAttribute(meta().id());
            return value == null ? defaultValue : valueOf(value);
        }
    }

    default T of(final VariantContext variant, final Genotype gt, final T defaultValue) {
        if (meta().target() != VCFAnnotationTarget.GENOTYPE) {
            throw new UnsupportedOperationException("wrong annotation target type");
        } else {
            final Object value = gt.getAnyAttribute(meta().id());
            return value == null ? defaultValue : valueOf(value);
        }
    }

    default T valueOf(final Object obj) {
        throw new UnsupportedOperationException("direct conversion without context is not supported");
    }


}
