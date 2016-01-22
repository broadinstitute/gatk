package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.api.services.genomics.model.VariantAnnotation;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;

import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;

/**
 * TODO: STUB
 */
public class VariantAnnotatorEngine {

    /**
     * Annotates the given variant context - adds all annotations that satisfy the predicate.
     */
    public VariantContext annotateContext(final FeatureContext features,
                                          final ReferenceContext ref,
                                          final VariantContext vc,
                                          final Map<String,PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap,
                                          final Predicate<VariantAnnotation> addAnnot) {
        return null;
    }

    public VariantContext annotateContext(final FeatureContext features,
                                          final ReferenceContext ref,
                                          final Map<String, AlignmentContext> stratifiedContexts,
                                          final VariantContext vc,
                                          final Map<String,PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        return null;
    }

    public VariantContext annotateContext(final FeatureContext features,
                                          final ReferenceContext ref,
                                          final Map<String, AlignmentContext> stratifiedContexts,
                                          final VariantContext vc) {
        return annotateContext(features, ref, stratifiedContexts, vc, null);
    }

    public Set<VCFHeaderLine> getVCFAnnotationDescriptions(){ return null;}

    public VariantContext annotateContextForActiveRegion(final ReferenceContext referenceContext,
                                                         final FeatureContext tracker,
                                                         final ReadLikelihoods<Allele> readAlleleLikelihoods,
                                                         final VariantContext call,
                                                         final boolean emitReferenceConfidence) {
        return null;
    }
}
