package org.broadinstitute.hellbender.utils.read.mergealignment;

/**
 * Given a set of alignments for a read or read pair, mark one alignment as primary, according to whatever
 * strategy is appropriate.  Any pre-existing primary designation is ignored, so if the aligner has selected an
 * appropriate primary alignment, this class should not be called.
 */
public interface PrimaryAlignmentSelectionStrategy {
    /**
     * When this method returns, one alignment has been marked as primary according to the implementation's strategy.
     *
     */
    void pickPrimaryAlignment(HitsForInsert hitsForInsert);
}
