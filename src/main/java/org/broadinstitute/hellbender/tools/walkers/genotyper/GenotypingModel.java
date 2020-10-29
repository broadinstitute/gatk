package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.dragstr.DragstrReferenceAnalyzer;

/**
 * A wrapping interface between the various versions of genotypers so as to keep them interchangeable.
 */
public interface GenotypingModel {
    @SuppressWarnings({"unchecked", "rawtypes"})
    public <A extends Allele> GenotypingLikelihoods<A> calculateLikelihoods(AlleleList<A> genotypingAlleles, GenotypingData<A> data, final byte[] paddedReference, final int offsetForRefIntoEvent, final DragstrReferenceAnalyzer dragstrs);

}
