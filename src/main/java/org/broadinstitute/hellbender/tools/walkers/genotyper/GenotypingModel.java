package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.dragstr.DragstrReferenceAnalyzer;

import java.util.Optional;

/**
 * A wrapping interface between the various versions of genotypers so as to keep them interchangeable.
 */
public interface GenotypingModel {
    @SuppressWarnings({"unchecked", "rawtypes"})
    public <A extends Allele> GenotypingLikelihoods<A> calculateLikelihoods(AlleleList<A> genotypingAlleles, GenotypingData<A> data,
        final byte[] paddedReference, final int offsetForRefIntoEvent, final DragstrReferenceAnalyzer dragstrs);

    default <A extends Allele> GenotypingLikelihoods<A> calculateLikelihoods(AlleleList<A> genotypingAlleles, GenotypingData<A> data,
        final byte[] paddedReference, final int offsetForRefIntoEvent, final DragstrReferenceAnalyzer dragstrs, Optional<GenotypingLikelihoods<A>> glsOverride) {
        if (glsOverride.isEmpty()) {
            return calculateLikelihoods(genotypingAlleles, data, paddedReference, offsetForRefIntoEvent, dragstrs);
        } else {
            throw new UnsupportedOperationException("This GenotypingModel class does not implement a GL-override mode.");
        }
    }
}
