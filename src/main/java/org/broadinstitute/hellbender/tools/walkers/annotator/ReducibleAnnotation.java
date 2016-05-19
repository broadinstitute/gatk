package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;

import java.util.Map;

/**
 * An interface for annotations that are calculated using raw data across samples, rather than the median (or median of median) of samples values
 */
public interface ReducibleAnnotation extends Annotation {
    public abstract String getRawKeyName();

    /**
     * Generate the raw data necessary to calculate the annotation. Raw data is the final endpoint for gVCFs.
     */
    public abstract Map<String, Object> annotateRawData(final ReferenceContext ref,
                                                        final VariantContext vc,
                                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap);

    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME
    public abstract void calculateRawData(VariantContext vc, Map<String, PerReadAlleleLikelihoodMap> pralm, ReducibleAnnotationData rawAnnotations);
}
