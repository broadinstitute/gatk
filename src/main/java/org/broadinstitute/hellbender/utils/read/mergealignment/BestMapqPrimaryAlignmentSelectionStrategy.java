package org.broadinstitute.hellbender.utils.read.mergealignment;

import htsjdk.samtools.SAMUtils;
import org.broadinstitute.hellbender.utils.read.mergealignment.HitsForInsert.NumPrimaryAlignmentState;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * This strategy was designed for TopHat output, but could be of general utility.  It picks the alignment with best MAPQ.
 * If paired-end, it is the alignment in which the sum of the MAPQs of both ends is the best.  In case of ties, one
 * is selected arbitrarily.  This strategy expects pair-aware alignments, with the corresponding alignment for each
 * mate of the pair correlated by HI (hit index) tag.  If the aligner has set a pair of alignments as primary, this
 * is used (assuming one of those alignments is not filtered out).  Otherwise the alignment pair with best MapQ is
 * selected.
 */
public final class BestMapqPrimaryAlignmentSelectionStrategy implements PrimaryAlignmentSelectionStrategy {
    // Give random number generator a seed so results are repeatable.  Used to pick a primary alignment from
    // multiple alignments with equal mapping quality.
    private final Random random = new Random(1);

    /**
     * Primary alignment was filtered out.  Need to select a new one.
     */
    @Override
    public void pickPrimaryAlignment(final HitsForInsert hits) {

        if (hits.numHits() == 0) throw new IllegalArgumentException("No alignments to pick from");
        hits.coordinateByHitIndex();
        // See if primary alignment is not already unambiguously determined.
        final NumPrimaryAlignmentState firstEndAlignmentState = hits.tallyPrimaryAlignments(true);
        final NumPrimaryAlignmentState secondEndAlignmentState = hits.tallyPrimaryAlignments(false);

        if ((firstEndAlignmentState == NumPrimaryAlignmentState.NONE && secondEndAlignmentState == NumPrimaryAlignmentState.NONE) ||
                firstEndAlignmentState == NumPrimaryAlignmentState.MORE_THAN_ONE ||
                secondEndAlignmentState == NumPrimaryAlignmentState.MORE_THAN_ONE) {
            // Need to use selected strategy for picking primary.

            // Find all the hits with the best MAPQ.
            final List<Integer> primaryAlignmentIndices = new ArrayList<>(hits.numHits());
            int bestMapQ = -1;
            for (int i = 0; i < hits.numHits(); ++i) {
                final int firstEndMapq;
                if (hits.getFirstOfPair(i) != null) {
                    firstEndMapq = hits.getFirstOfPair(i).getMappingQuality();
                } else {
                    firstEndMapq = 0;
                }
                final int secondEndMapq;
                if (hits.getSecondOfPair(i) != null) {
                    secondEndMapq = hits.getSecondOfPair(i).getMappingQuality();
                } else {
                    secondEndMapq = 0;
                }
                int thisMapQ = SAMUtils.combineMapqs(firstEndMapq, secondEndMapq);
                if (thisMapQ > bestMapQ) {
                    bestMapQ = thisMapQ;
                    primaryAlignmentIndices.clear();
                }
                if (thisMapQ == bestMapQ) primaryAlignmentIndices.add(i);
            }

            // Of all the hits with the best MAPQ, randomly select one to be primary.
            final int primaryAlignmentIndex;
            if (primaryAlignmentIndices.size() == 1) primaryAlignmentIndex = primaryAlignmentIndices.get(0);
            else if (primaryAlignmentIndices.size() > 1) primaryAlignmentIndex =
                    primaryAlignmentIndices.get(random.nextInt(primaryAlignmentIndices.size()));
            else throw new IllegalStateException("Never found a best MAPQ -- should never happen");

            hits.setPrimaryAlignment(primaryAlignmentIndex);
        }
    }
}
