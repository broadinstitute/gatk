package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Map;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.CPX_SV_SYB_ALT_ALLELE_STR;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR;

final class CpxVariantType extends SvType {

    @Override
    public final String toString() {
        return CPX_SV_SYB_ALT_ALLELE_STR;
    }

    CpxVariantType(final SimpleInterval affectedRefRegion, final int altHaplotypeSequenceLength,
                   final Map<String, String> typeSpecificExtraAttributes) {
        super(getIDString(affectedRefRegion),
                Allele.create(SimpleSVType.createBracketedSymbAlleleString(CPX_SV_SYB_ALT_ALLELE_STR)),
                altHaplotypeSequenceLength - affectedRefRegion.size(), typeSpecificExtraAttributes);
    }

    private static String getIDString(final SimpleInterval affectedRefRegion) {
        return CPX_SV_SYB_ALT_ALLELE_STR + INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                affectedRefRegion.toString();
    }
}
