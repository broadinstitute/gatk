package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.samtools.util.Locatable;
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
        final byte[] paddedReference, final int offsetForRefIntoEvent, final DragstrReferenceAnalyzer dragstrs, final Locatable eventLocus);

    default <A extends Allele> GenotypingLikelihoods<A> calculateLikelihoods(AlleleList<A> genotypingAlleles, GenotypingData<A> data,
        final byte[] paddedReference, final int offsetForRefIntoEvent, final DragstrReferenceAnalyzer dragstrs,
        Optional<GenotypingLikelihoods<A>> genotypeLikelihoodsOverride, Optional<GenotypingLikelihoods<A>> genotypePosteriorsOverride,
                                                                             final Locatable eventLocus) {
        if (genotypeLikelihoodsOverride.isEmpty()) {
            return calculateLikelihoods(genotypingAlleles, data, paddedReference, offsetForRefIntoEvent, dragstrs, eventLocus);
        } else {
            throw new UnsupportedOperationException("This GenotypingModel class does not implement a GL-override mode.");
        }
    }
}
