package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

enum VariantTypeMode {
    SNP,
    INDEL,
    BOTH;

    static boolean checkVariationClass(final VariantContext vc,
                                       final VariantContext resourceVC) {
        switch (resourceVC.getType()) {
            case SNP:
            case MNP:
                return checkVariationClass(vc, SNP);
            case INDEL:
            case MIXED:
            case SYMBOLIC:
                return checkVariationClass(vc, INDEL);
            default:
                return false;
        }
    }

    static boolean checkVariationClass(final VariantContext vc,
                                       final Allele allele,
                                       final VariantTypeMode mode) {
        switch (mode) {
            case SNP:
                //note that spanning deletions are considered SNPs by this logic
                return vc.getReference().length() == allele.length();
            case INDEL:
                return (vc.getReference().length() != allele.length()) || allele.isSymbolic();
            case BOTH:
                return true;
            default:
                throw new IllegalStateException("Encountered unknown mode: " + mode);
        }
    }

    static boolean checkVariationClass(final VariantContext vc,
                                       final VariantTypeMode mode) {
        switch (mode) {
            case SNP:
                return vc.isSNP() || vc.isMNP();
            case INDEL:
                return vc.isStructuralIndel() || vc.isIndel() || vc.isMixed() || vc.isSymbolic();
            case BOTH:
                return true;
            default:
                throw new IllegalStateException("Encountered unknown mode: " + mode);
        }
    }
}
