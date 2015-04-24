package org.broadinstitute.hellbender.tools.picard.vcf.filter;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * Filters out a record if all variant samples have depth lower than the given value.
 */
public final class DepthFilter implements GenotypeFilter {
    private final int minDepth;

    public DepthFilter(final int minDepth) {
        this.minDepth = minDepth;
    }

    @Override
    public String filter(final VariantContext ctx, final Genotype gt) {
        if (gt.getDP() < minDepth) return "LowDP";
        else return null;
   }
}
