package org.broadinstitute.hellbender.tools.walkers.annotator;

public interface VCFAnnotationReducer {

    void reduceAll(final VCFAnnotationContext target, final VCFAnnotationContext ... input);

    <T> void reduce(final VCFAnnotation<T> annotation, final VCFAnnotationContext target, final VCFAnnotationContext ... input);
}
