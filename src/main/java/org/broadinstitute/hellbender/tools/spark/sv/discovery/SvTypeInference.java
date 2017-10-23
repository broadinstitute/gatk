package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.exceptions.GATKException;

class SvTypeInference {

    /**
     * @return inferred type of variant (as listed in {@link SvType}) based on input {@link NovelAdjacencyReferenceLocations}.
     */
    @VisibleForTesting
    static SvType inferFromNovelAdjacency(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations) {

        final int start = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd();
        final int end = novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart();
        final StrandSwitch strandSwitch = novelAdjacencyReferenceLocations.strandSwitch;

        final SvType type;
        if (strandSwitch == StrandSwitch.NO_SWITCH) { // no strand switch happening, so no inversion

            final boolean hasDupAnnotation = novelAdjacencyReferenceLocations.complication.hasDuplicationAnnotation();
            final boolean hasNoInsertedSeq = novelAdjacencyReferenceLocations.complication.getInsertedSequenceForwardStrandRep().isEmpty();

            if (hasDupAnnotation) { // could be expansion or contraction
                final boolean isExpansion = novelAdjacencyReferenceLocations.complication.getDupSeqRepeatNumOnCtg()
                                            > novelAdjacencyReferenceLocations.complication.getDupSeqRepeatNumOnRef();
                if (isExpansion) {
                    // simple expansion of repeat 1 -> 2 with possible inserted seq in between, or complex expansion
                    type = new SimpleSVType.DuplicationTandem(novelAdjacencyReferenceLocations);
                } else {
                    // simple contraction of repeat 2 -> 1, or complex contraction
                    type = new SimpleSVType.Deletion(novelAdjacencyReferenceLocations);
                }
            } else { // simple ins del
                if (start == end) { // must be for simple ins records
                    if (hasNoInsertedSeq) {
                        throw new GATKException("Something went wrong in type inference, there's suspected insertion happening but no inserted sequence could be inferred "
                                + novelAdjacencyReferenceLocations.toString());
                    } else {
                        type = new SimpleSVType.Insertion(novelAdjacencyReferenceLocations); // simple insertion
                    }
                } else {
                    // TODO: 9/27/17 wait for the VCF spec to iron out how to distinguish the two
                    if (hasNoInsertedSeq) {
                        type = new SimpleSVType.Deletion(novelAdjacencyReferenceLocations); // clean deletion
                    } else {
                        type = new SimpleSVType.Deletion(novelAdjacencyReferenceLocations); // scarred deletion
                    }
                }
            }
        } else {
            type = new SimpleSVType.Inversion(novelAdjacencyReferenceLocations);
        }

        // developer check to make sure new types are treated correctly
        try {
            SimpleSVType.TYPES.valueOf(type.toString());
        } catch (final IllegalArgumentException ex) {
            throw new GATKException.ShouldNeverReachHereException("Inferred type is not known yet: " + type.toString(), ex);
        }

        return type;
    }
}
