package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;

import java.util.List;

public final class SVMappingQualityFilter implements StructuralVariantFilter {
    static final long serialVersionUID = 1L;

    private static final String attributeKey = GATKSVVCFConstants.MAPPING_QUALITIES;
    private final int threshold;

    public SVMappingQualityFilter( final int threshold) {
        this.threshold = threshold;
    }

    @Override
    public String getName() {
        return GATKSVVCFConstants.ASSEMBLY_BASED_VARIANT_MQ_FILTER_KEY;
    }

    @Override
    public boolean test(final VariantContext variantContext) {
        if ( !variantContext.hasAttribute(GATKSVVCFConstants.CONTIG_NAMES) )
            return true;

        final List<String> mapQuals = SVUtils.getAttributeAsStringList(variantContext, attributeKey);
        int maxMQ = 0;
        for (final String mapQual : mapQuals) {
            Integer integer = Integer.valueOf(mapQual);
            maxMQ = maxMQ < integer ? integer : maxMQ;
        }
        return maxMQ >= threshold;
    }
}
