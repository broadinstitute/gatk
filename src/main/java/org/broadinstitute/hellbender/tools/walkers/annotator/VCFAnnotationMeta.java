package org.broadinstitute.hellbender.tools.walkers.annotator;

public @interface VCFAnnotationMeta {

    VCFAnnotationTarget target();

    String id();

}
