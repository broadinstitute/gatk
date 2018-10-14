package org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation;

import java.util.Arrays;
import java.util.List;

public final class CopyNumberAnnotations {
    public static AnnotationKey<Double> GC_CONTENT = new AnnotationKey<>(
            "GC_CONTENT",
            Double.class,
            gcContent -> (0. <= gcContent && gcContent <= 1.) || Double.isNaN(gcContent));

    public static AnnotationKey<Double> MAPPABILITY = new AnnotationKey<>(
            "MAPPABILITY",
            Double.class,
            mappability -> (0. <= mappability && mappability <= 1.) || Double.isNaN(mappability));

    public static AnnotationKey<Double> SEGMENTAL_DUPLICATION_CONTENT = new AnnotationKey<>(
            "SEGMENTAL_DUPLICATION_CONTENT",
            Double.class,
            segmentalDuplicationContent -> (0. <= segmentalDuplicationContent && segmentalDuplicationContent <= 1.) || Double.isNaN(segmentalDuplicationContent));

    /**
     * This defines the canonical order of these annotations.
     */
    public static List<AnnotationKey<?>> ANNOTATIONS = Arrays.asList(GC_CONTENT, MAPPABILITY, SEGMENTAL_DUPLICATION_CONTENT);
}
