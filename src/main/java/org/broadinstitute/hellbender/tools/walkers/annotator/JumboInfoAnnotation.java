package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * INFO annotations that look at more inputs than regular annotations
 */
public interface JumboInfoAnnotation extends VariantAnnotation{

    default VCFCompoundHeaderLine.SupportedHeaderLineType annotationType() { return VCFCompoundHeaderLine.SupportedHeaderLineType.INFO; }

    Map<String, Object> annotate(final ReferenceContext ref,
                                                 final FeatureContext features,
                                                 final VariantContext vc,
                                                 final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                                                 final AlleleLikelihoods<Fragment, Allele> fragmentLikelihoods,
                                                 final AlleleLikelihoods<Fragment, Haplotype> haplotypeLikelihoods);
}