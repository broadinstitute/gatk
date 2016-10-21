package org.broadinstitute.hellbender.utils.read.mergealignment;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CoordMath;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * When it is necessary to pick a primary alignment from a group of alignments for a read, pick the one that maps
 * the earliest base in the read. This implementation only works for fragments, not for pairs.
 * If there are multiple alignments that all start mapping at the same offest in the read, pick the one with the best
 * MAPQ.  If there are multiple alignments that have the earliest mapping and that have the same MAPQ, pick one randomly.
 */
public final class EarliestFragmentPrimaryAlignmentSelectionStrategy implements PrimaryAlignmentSelectionStrategy {
    // Give random number generator a seed so results are repeatable.  Used to pick a primary alignment from
    // multiple alignments with equal mapping quality.
    private final Random random = new Random(1);

    @Override
    public void pickPrimaryAlignment(final HitsForInsert hitsForInsert) {

        if (hitsForInsert.numHits() == 0) throw new IllegalArgumentException("No alignments to pick from");

        // Gather the earliest alignment(s) with best MAPQ
        final List<Integer> earliestAlignments = new ArrayList<>();
        int earliestMappedBase = Integer.MAX_VALUE;
        int bestMapQ = -1;
        for (int i = 0; i < hitsForInsert.numHits(); ++i) {
            final SAMRecord rec = hitsForInsert.getFragment(i);
            if (rec.getReadUnmappedFlag()) continue;
            final int thisFirstMappedBase = getIndexOfFirstAlignedBase(rec);
            final int thisMapQ = rec.getMappingQuality();
            if (thisFirstMappedBase < earliestMappedBase ||
                    (thisFirstMappedBase == earliestMappedBase && thisMapQ > bestMapQ)) {
                earliestAlignments.clear();
                earliestAlignments.add(i);
                earliestMappedBase = thisFirstMappedBase;
                bestMapQ = thisMapQ;
            } else if (thisFirstMappedBase == earliestMappedBase && thisMapQ == bestMapQ) {
                earliestAlignments.add(i);
            } // else it is not the earliest or the best so skip it.
        }


        if (earliestAlignments.size() == 1) {
            // If only one best, pick it.
            hitsForInsert.setPrimaryAlignment(earliestAlignments.get(0));
        } else {
            // Arbitrarily select one of the best
            hitsForInsert.setPrimaryAlignment(earliestAlignments.get(random.nextInt(earliestAlignments.size())));
        }
    }

    /**
     * Returns 1-based index of first base in read that corresponds to M in CIGAR string.
     * Note that first is relative to 5' end, so that for reverse-strand alignment, the index of
     * the last base aligned is computed relative to the end of the read.
     */
    int getIndexOfFirstAlignedBase(final SAMRecord rec) {
        final List<AlignmentBlock> alignmentBlocks = rec.getAlignmentBlocks();
        if (rec.getReadNegativeStrandFlag()) {
            final AlignmentBlock alignmentBlock = alignmentBlocks.get(alignmentBlocks.size() - 1);
            return rec.getReadLength() - CoordMath.getEnd(alignmentBlock.getReadStart(), alignmentBlock.getLength()) + 1;
        } else {
            return alignmentBlocks.get(0).getReadStart();
        }
    }
}
