package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Map;

/**
 * Annotations relevant to the INFO field of the variant file (ie annotations for sites).
 */
public interface InfoFieldAnnotation extends VariantAnnotation {

    default AnnotationType annotationType() { return AnnotationType.INFO; }

    /**
     * Computes the annotation for the given variant and the likelihoods per read.
     * Returns a map from annotation keys to values (may be empty if no annotation is to be added).
     *
     * @param ref Reference context, may be null
     * @param vc Variant to be annotated. Not null.
     * @param likelihoods likelihoods indexed by sample, allele, and read within sample
     */
    public Map<String, Object> annotate(final ReferenceContext ref,
                                                 final VariantContext vc,
                                                 final AlleleLikelihoods<GATKRead, Allele> likelihoods);
}