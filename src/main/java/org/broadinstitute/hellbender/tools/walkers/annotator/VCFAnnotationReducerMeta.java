package org.broadinstitute.hellbender.tools.walkers.annotator;

public @interface VCFAnnotationReducerMeta {
    String[] value() default {};
    String[] requires() default {};
}
