package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;

public interface GenotypersModel {
    @SuppressWarnings({"unchecked", "rawtypes"})
    <A extends Allele> GenotypingLikelihoods<A> calculateLikelihoods(AlleleList<A> genotypingAlleles, GenotypingData<A> data);
}
