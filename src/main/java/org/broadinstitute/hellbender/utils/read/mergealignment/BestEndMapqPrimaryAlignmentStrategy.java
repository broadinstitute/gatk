package org.broadinstitute.hellbender.utils.read.mergealignment;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;

import java.io.Serializable;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

/**
 * For an aligner that aligns each end independently, select the alignment for each end with the best MAPQ, and
 * make that the primary.  The primary alignments are then correlated so that their mate info points to each
 * other, but all non-primary alignments are uncorrelated.
 */
public final class BestEndMapqPrimaryAlignmentStrategy implements PrimaryAlignmentSelectionStrategy {
    private static final MapqComparator MAPQ_COMPARATOR = new MapqComparator();

    // Give random number generator a seed so results are repeatable.  Used to pick a primary alignment from
    // multiple alignments with equal mapping quality.
    private final Random random = new Random(1);

    /**
     * Primary alignment was filtered out.  Need to select a new one.
     */
    @Override
    public void pickPrimaryAlignment(final HitsForInsert hits) {

        if (hits.numHits() == 0) throw new IllegalArgumentException("No alignments to pick from");
        Collections.sort(hits.firstOfPairOrFragment, MAPQ_COMPARATOR);
        Collections.sort(hits.secondOfPair, MAPQ_COMPARATOR);

        randomlySelectPrimaryFromBest(hits.firstOfPairOrFragment);
        randomlySelectPrimaryFromBest(hits.secondOfPair);
        hits.setPrimaryAlignment(0);

        if (!hits.isPaired()) return;

        // For non-primary alignments, de-correlate them so that the mate fields don't point at some
        // arbitrary alignment for the other end.

        // No non-primary alignments for one of the ends so nothing to do.
        if (hits.firstOfPairOrFragment.size() <= 1 || hits.secondOfPair.size() <= 1) return;
        final int amountToSlide = hits.firstOfPairOrFragment.size() - 1;
        for (int i = 0; i < amountToSlide; ++i) {
            hits.secondOfPair.add(1, null);
        }
    }


    /**
     * Randomly picks one of the best alignments and puts it into the 0th slot of the list.
     * @param recs List of alignments sorted in descending order of alignment quality.
     */
    private void randomlySelectPrimaryFromBest(List<SAMRecord> recs) {
        if (recs.isEmpty()) return;
        final int bestMapq = recs.get(0).getMappingQuality();
        int i;
        for (i = 1; i < recs.size() && recs.get(i).getMappingQuality() == bestMapq; ++i) {
        }
        final int bestIndex = random.nextInt(i);
        if (bestIndex == 0) return;
        final SAMRecord tmp = recs.get(0);
        recs.set(0, recs.get(bestIndex));
        recs.set(bestIndex, tmp);
    }

    // Sorts in descending order, but 255 is considered > 0 but < 1, and unmapped is worst of all
    private static class MapqComparator implements Comparator<SAMRecord>, Serializable {
        private static final long serialVersionUID = 6763153425070516820L;

        @Override
        public int compare(final SAMRecord rec1, final SAMRecord rec2) {
            if (rec1.getReadUnmappedFlag()) {
                if (rec2.getReadUnmappedFlag()) return 0;
                else return 1;
            } else if (rec2.getReadUnmappedFlag()) {
                return -1;
            }
            return -SAMUtils.compareMapqs(rec1.getMappingQuality(), rec2.getMappingQuality());
        }
    }
}
