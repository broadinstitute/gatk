package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * FORMAT annotations that look at more inputs than regular annotations
 */
public interface JumboGenotypeAnnotation extends VariantAnnotation{

    default AnnotationType annotationType() { return AnnotationType.GENOTYPE; }

    void annotate(final ReferenceContext ref,
                                  final FeatureContext features,
                                  final VariantContext vc,
                                  final Genotype g,
                                  final GenotypeBuilder gb,
                                  final AlleleLikelihoods<GATKRead, Allele> readLikelihoods,
                                  final AlleleLikelihoods<Fragment, Allele> fragmentLikelihoods,
                                  final AlleleLikelihoods<Fragment, Haplotype> haplotypeLikelihoods);
}