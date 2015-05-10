package org.broadinstitute.hellbender.tools.picard.vcf.filter;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * Genotype filter that filters out genotypes below a given quality threshold.
 *
 * @author tfennell
 */
public final class GenotypeQualityFilter implements GenotypeFilter {
    private final int minGq;

    public GenotypeQualityFilter(final int minGq) {
        this.minGq = minGq;
    }

    @Override
    public String filter(final VariantContext ctx, final Genotype gt) {
        if (gt.getGQ() < minGq) return "LowGQ";
        else return null;
    }
}
