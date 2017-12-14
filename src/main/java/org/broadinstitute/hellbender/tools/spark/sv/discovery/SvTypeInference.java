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
            if (start==end) { // something is inserted
                final boolean hasNoDupSeq = !novelAdjacencyReferenceLocations.complication.hasDuplicationAnnotation();
                final boolean hasNoInsertedSeq = novelAdjacencyReferenceLocations.complication.getInsertedSequenceForwardStrandRep().isEmpty();
                if (hasNoDupSeq) {
                    if (hasNoInsertedSeq) {
                        throw new GATKException("Something went wrong in type inference, there's suspected insertion happening but no inserted sequence could be inferred "
                                + novelAdjacencyReferenceLocations.toString());
                    } else {
                        type = new SimpleSVType.Insertion(novelAdjacencyReferenceLocations); // simple insertion (no duplication)
                    }
                } else {
                    if (hasNoInsertedSeq) {
                        type = new SimpleSVType.DuplicationTandem(novelAdjacencyReferenceLocations); // clean expansion of repeat 1 -> 2, or complex expansion
                    } else {
                        type = new SimpleSVType.DuplicationTandem(novelAdjacencyReferenceLocations); // expansion of 1 repeat on ref to 2 repeats on alt with inserted sequence in between the 2 repeats
                    }
                }
            } else {
                final boolean hasNoDupSeq = !novelAdjacencyReferenceLocations.complication.hasDuplicationAnnotation();
                final boolean hasNoInsertedSeq = novelAdjacencyReferenceLocations.complication.getInsertedSequenceForwardStrandRep().isEmpty();
                if (hasNoDupSeq) {
                    if (hasNoInsertedSeq) {
                        type = new SimpleSVType.Deletion(novelAdjacencyReferenceLocations); // clean deletion
                    } else {
                        type = new SimpleSVType.Deletion(novelAdjacencyReferenceLocations); // scarred deletion
                    }
                } else {
                    if (hasNoInsertedSeq) {
                        type = new SimpleSVType.Deletion(novelAdjacencyReferenceLocations); // clean contraction of repeat 2 -> 1, or complex contraction
                    } else {
                        throw new GATKException("Something went wrong in type inference, there's suspected deletion happening but both inserted sequence and duplication exits (not supported yet): "
                                + novelAdjacencyReferenceLocations.toString());
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
