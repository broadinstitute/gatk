package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;

import java.util.Map;

/**
 * An interface for annotations that are calculated using raw data across samples, rather than the median (or median of median) of samples values
 *
 * TODO: this interface will eventually contain the actual 'reduce' functionality.
 */
public interface ReducibleAnnotation extends Annotation {
    public abstract String getRawKeyName();

    /**
     * Generate the raw data necessary to calculate the annotation. Raw data is the final endpoint for gVCFs.
     * @param ref the reference context for this annotation
     * @param vc the variant context to annotate
     * @param stratifiedPerReadAlleleLikelihoodMap map of: sampleName -> read -> allele -> likelihood or null if there's no reads to process.
     */
    public abstract Map<String, Object> annotateRawData(final ReferenceContext ref,
                                                        final VariantContext vc,
                                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap);

    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME: there's no way to properly use generics here because of method overriding
    public abstract void calculateRawData(VariantContext vc, Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap, ReducibleAnnotationData rawAnnotations);
}
