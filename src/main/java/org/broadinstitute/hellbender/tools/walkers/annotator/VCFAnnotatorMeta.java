package org.broadinstitute.hellbender.tools.walkers.annotator;

public @interface VCFAnnotatorMeta {
    String[] value() default {};
    String[] requiredGenotypeAnnotations() default {};
    String[] requiredVariantAnnotations() default {};
    boolean requiresGenotypes() default false;
    boolean requiresAlignment() default false;
    boolean requiresReference() default false;
    boolean requiresLikelihoods() default false;
    int requiredReferenceWindow() default 0;
    double variance() default 0; // exact.
    double bias() default 0; // no bias.

}
