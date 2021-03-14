package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;

public final class SVAlignmentLengthFilter implements StructuralVariantFilter {
    static final long serialVersionUID = 1L;

    private static final String attributeKey = GATKSVVCFConstants.MAX_ALIGN_LENGTH;
    private final int threshold;

    public SVAlignmentLengthFilter( final int threshold) {
        this.threshold = threshold;
    }

    @Override
    public String getName() {
        return GATKSVVCFConstants.ASSEMBLY_BASED_VARIANT_ALN_LENGTH_FILTER_KEY;
    }

    @Override
    public boolean test(final VariantContext variantContext) {
        if ( !variantContext.hasAttribute(GATKSVVCFConstants.CONTIG_NAMES) )
            return true;

        final int alnLength = variantContext.getAttributeAsInt(attributeKey, 0);
        return alnLength >= threshold;
    }
}
