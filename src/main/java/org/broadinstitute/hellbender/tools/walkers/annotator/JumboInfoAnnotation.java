package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Map;

/**
 * INFO annotations that look at more inputs than regular annotations
 */
public interface JumboInfoAnnotation extends VariantAnnotation{

    default AnnotationType annotationType() { return AnnotationType.INFO; }

    Map<String, Object> annotate(final ReferenceContext ref,
                                                 final FeatureContext features,
                                                 final VariantContext vc,
                                                 final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                                                 final AlleleLikelihoods<Fragment, Allele> fragmentLikelihoods,
                                                 final AlleleLikelihoods<Fragment, Haplotype> haplotypeLikelihoods);
}