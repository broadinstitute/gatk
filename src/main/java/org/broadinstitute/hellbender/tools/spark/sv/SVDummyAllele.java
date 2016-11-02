package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * Simply a dummy class so that {@link org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods},
 * which takes an {@link Allele} type parameter, won't complain.
 */
class SVDummyAllele extends Allele{
    private static final long serialVersionUID = 1L;

    static final SVDummyAllele SymbolicReferenceAllele = new SVDummyAllele("N".getBytes(), true);
    static final SVDummyAllele SymbolicAlternateAllele = new SVDummyAllele("A".getBytes(), false);

    SVDummyAllele(final VariantContext vc, final boolean isRef){
        super(vc.getAlleles().get(0).getBases(), isRef);
    }

    SVDummyAllele(final byte[] bases, final boolean isRef){
        super(bases, isRef);
    }

    SVDummyAllele(final Allele allele) { super(allele, false);}
}
