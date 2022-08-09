package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * This code and logic for determining variant types was mostly retained from VQSR.
 * Note that there may be some inconsistencies and room for improvement in these definitions;
 * see comments in https://github.com/broadinstitute/gatk/pull/7954.
 */
public enum VariantType {
    SNP,
    INDEL;

    /**
     * Returns true if both {@code vc} and {@code resourceVC} are the same variant type,
     * following our definitions.
     */
    public static boolean checkVariantType(final VariantContext vc,
                                           final VariantContext resourceVC) {
        switch (resourceVC.getType()) {
            case SNP:
            case MNP:
                return getVariantType(vc) == SNP;
            case INDEL:
            case MIXED:
            case SYMBOLIC:
                return getVariantType(vc) == INDEL;
            default:
                return false;
        }
    }

    public static VariantType getVariantType(final VariantContext vc) {
        if (vc.isSNP() || vc.isMNP()) {
            return SNP;
        } else if (vc.isStructuralIndel() || vc.isIndel() || vc.isMixed() || vc.isSymbolic()) {
            return INDEL;
        } else {
            throw new IllegalStateException("Encountered unknown variant type: " + vc.getType());
        }
    }

    /**
     * Note that spanning deletions are expected to be filtered out upstream of this method
     * to preserve VQSR behavior; we do not explicitly check this.
     * See VariantDataManager#checkVariationClass(VariantContext, Allele, VariantRecalibratorArgumentCollection.Mode),
     * from which this method originated.
     */
    public static VariantType getAlleleSpecificVariantType(final VariantContext vc,
                                                           final Allele allele) {
        if (vc.getReference().length() == allele.length()) {
            // note that spanning deletions would be considered SNPs by this logic
            return SNP;
        }
        return INDEL;
    }
}
