package org.broadinstitute.hellbender.tools.picard.vcf.filter;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * An interface for classes that perform Genotype filtration. Implementations are expected to take in a VariantContext
 * and a single Genotype and return either null (for no filter) or a specific filter string.
 *
 * @author Tim Fennell
 */
public interface GenotypeFilter {
    /** Test whether or not the genotype should be filtered out. If so return a filter string, otherwise return null. */
    public String filter(final VariantContext ctx, final Genotype gt);
}
