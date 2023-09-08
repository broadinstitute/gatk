package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
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

    default VCFCompoundHeaderLine.SupportedHeaderLineType annotationType() { return VCFCompoundHeaderLine.SupportedHeaderLineType.INFO; }

    // TODO: fix ? extends Locatable - that actually is not needed, just pass the allele parts of the likelihoods
    Map<String, Object> annotate(final ReferenceContext ref,
                                 final FeatureContext features,
                                 final VariantContext vc,
                                 final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                                 final AlleleLikelihoods<? extends Locatable, Allele> fragmentLikelihoods,
                                 final AlleleLikelihoods<? extends Locatable, Haplotype> haplotypeLikelihoods);
}