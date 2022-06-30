package org.broadinstitute.hellbender.tools.sv.concordance;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVLinkage;

public class SVConcordanceLinkage extends CanonicalSVLinkage<SVCallRecord> {

    public SVConcordanceLinkage(final SAMSequenceDictionary dictionary) {
        super(dictionary, false);
    }

    @Override
    public boolean areClusterable(final SVCallRecord a, final SVCallRecord b) {
        final StructuralVariantType aType = a.getType();
        final StructuralVariantType bType = b.getType();
        // Don't allow CNV/DEL or CNV/DUP matching, which is problematic for concordance calculations
        if ((aType == StructuralVariantType.CNV || bType != StructuralVariantType.CNV) && aType != bType) {
            return false;
        }
        return super.areClusterable(a, b);
    }
}
