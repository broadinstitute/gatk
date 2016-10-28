package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;

import java.util.List;

/**
 * Represents an annotation that is computed for a single genotype.
 */
public abstract class GenotypeAnnotation extends VariantAnnotation{

    /**
     * Computes the annotation for the given genotype and the likelihoods per read.
     * Expected to modified the passed genotype builder.
     *
     * @param ref Reference context, may be null
     * @param vc Variant to be annotated. Not null.
     * @param likelihoods matrix of likelihoods indexed by allele and read
     * @param g the genotype to annotate. May be null.
     * @param gb the builder to modify and annotations to. Not null.
     */
    public abstract void annotate(final ReferenceContext ref,
                                  final VariantContext vc,
                                  final Genotype g,
                                  final GenotypeBuilder gb,
                                  final ReadLikelihoods<Allele> likelihoods);

    /**
     * Return the descriptions used for the VCF FORMAT meta field.
     * Subclasses must ensure that this list is not null and does not contain null.
     */
    public abstract List<VCFFormatHeaderLine> getDescriptions();
}