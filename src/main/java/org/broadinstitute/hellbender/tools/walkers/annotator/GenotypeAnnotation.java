package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableList;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Represents an annotation that is computed for a single genotype.
 */
public interface GenotypeAnnotation extends VariantAnnotation {

    default VCFCompoundHeaderLine.SupportedHeaderLineType annotationType() { return VCFCompoundHeaderLine.SupportedHeaderLineType.FORMAT; }

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
    public void annotate(final ReferenceContext ref,
                                  final VariantContext vc,
                                  final Genotype g,
                                  final GenotypeBuilder gb,
                                  final AlleleLikelihoods<GATKRead, Allele> likelihoods);
}