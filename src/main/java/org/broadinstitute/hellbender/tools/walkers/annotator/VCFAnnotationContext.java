package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;

public interface VCFAnnotationContext {

    VCFAnnotationDictionary dictionary();

    SampleList samples();

    int position();

    String chromosome();

    Allele reference();

    VariantContext toVariant();


    <T> T get(final VCFAnnotation<? extends T> annotation);

    <T> void set(final VCFAnnotation<? extends T> annotation, T value);

    boolean has(final VCFAnnotation<?> annotation);

    default <T> T get(final String key, final Class<T> clazz) {
        final VCFAnnotation<T> annotation = dictionary().forVariant(key, clazz);
        return annotation == null ? null : get(annotation);
    }

    <T> T get(final VCFAnnotation<? extends T> annotation, final String sample, final Class<T> clazz);


    void clear();

    void remove(final VCFAnnotation<?> annotation);

    AlignmentContext alignment();

    default boolean hasAlignment() {
        return alignment() != null;
    }

    ReferenceContext reference();

    default boolean hasReference() {
        return reference() != null;
    }

    AlleleList<?> alleles();

    default boolean satisfiesRequirements(final VCFAnnotatorMeta annotatorMeta) {
        if (annotatorMeta.requiresAlignment() && !hasAlignment()) {
            return false;
        } else if (annotatorMeta.requiresGenotypes() && !hasGenotypes()) {
            return false;
        } else if (annotatorMeta.requiresReference() && !hasReference()) {
            return false;
        } else if (annotatorMeta.requiresLikelihoods() && !hasLikelihoods()) {
            return false;
        } else {
            for (final String key : annotatorMeta.requiredVariantAnnotations()) {

            }
        }
     }
}
