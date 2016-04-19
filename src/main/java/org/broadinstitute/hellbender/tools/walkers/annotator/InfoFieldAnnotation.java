package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
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
     * Returns a map from annotation keys to values.
     *
     * @param ref Reference context, may be null
     * @param vc Variant to be annotated. Not null.
     * @param stratifiedPerReadAlleleLikelihoodMap map of likelihoods per read. May be null.
     */
    public abstract Map<String, Object> annotate(final ReferenceContext ref,
                                                 final VariantContext vc,
                                                 final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap);

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