package org.broadinstitute.hellbender.tools.walkers.annotator;

public interface VCFAnnotationDictionary {


    <T> VCFAnnotation<T> forVariant(String key, Class<T> clazz);

    public VCFAnnotation
}
