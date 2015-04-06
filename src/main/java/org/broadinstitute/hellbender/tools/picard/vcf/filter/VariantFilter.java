package org.broadinstitute.hellbender.tools.picard.vcf.filter;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFilterHeaderLine;

import java.util.*;

/**
 * Interface for classes that can generate filters for VariantContexts. The contract is that a
 * VariantContext is provided, and if the variant should be filtered out then the filter string
 * should be returned, otherwise null.
 *
 * @author tfennell
 */
public interface VariantFilter {
    /** Check to see if the VariantContext should have a filter applied to it. If so return the filter string, otherwise return null. */
    public String filter(final VariantContext ctx);

    /** Return VCF header lines that define filters that may be applied by the VariantFilter. */
    public List<VCFFilterHeaderLine> headerLines();
}
