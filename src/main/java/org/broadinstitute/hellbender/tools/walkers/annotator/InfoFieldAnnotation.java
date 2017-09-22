package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotation;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Annotations relevant to the INFO field of the variant file (ie annotations for sites).
 */
public abstract class InfoFieldAnnotation extends VariantAnnotation{

    /**
     * Computes the annotation for the given variant and the likelihoods per read.
     * Returns a map from annotation keys to values (may be empty if no annotation is to be added).
     *
     * @param ref Reference context, may be null
     * @param vc Variant to be annotated. Not null.
     * @param likelihoods likelihoods indexed by sample, allele, and read within sample
     */
    public abstract Map<String, Object> annotate(final ReferenceContext ref,
                                                 final VariantContext vc,
                                                 final ReadLikelihoods<Allele> likelihoods);

    /**
     * Returns the descriptions used for the VCF INFO meta field.
     * Subclasses must ensure that this list is not null and does not contain null.
     */
    public List<VCFInfoHeaderLine> getDescriptions() {
        final List<VCFInfoHeaderLine> lines = new ArrayList<>(getKeyNames().size());
        for (final String key : getKeyNames()) {
            lines.add(GATKVCFHeaderLines.getInfoLine(key));
        }
        return lines;
    }
}