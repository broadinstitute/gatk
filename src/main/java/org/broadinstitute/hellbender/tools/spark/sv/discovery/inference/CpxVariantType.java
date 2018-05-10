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

    CpxVariantType(final SimpleInterval affectedRefRegion, final byte[] refBases, final int altHaplotypeSequenceLength,
                   final Map<String, Object> typeSpecificExtraAttributes) {
        super(affectedRefRegion.getContig(),
                affectedRefRegion.getStart(),
                affectedRefRegion.getEnd(),
                getIDString(affectedRefRegion), //id
                Allele.create(refBases, true),
                Allele.create(SimpleSVType.createBracketedSymbAlleleString(CPX_SV_SYB_ALT_ALLELE_STR)), // alt allele
                altHaplotypeSequenceLength - affectedRefRegion.size(), // svlen
                typeSpecificExtraAttributes);
    }

    @Override
    public final boolean hasApplicableEnd() {
        return true;
    }
    @Override
    public final boolean hasApplicableLength() {
        return true;
    }

    private static String getIDString(final SimpleInterval affectedRefRegion) {
        return CPX_SV_SYB_ALT_ALLELE_STR + INTERVAL_VARIANT_ID_FIELD_SEPARATOR +
                affectedRefRegion.toString();
    }
}
