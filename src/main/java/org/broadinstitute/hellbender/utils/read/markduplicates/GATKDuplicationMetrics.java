package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import picard.sam.DuplicationMetrics;

import java.io.Serializable;

/**
 * Metrics that are calculated during the process of marking duplicates
 * within a stream of SAMRecords.
 */
public final class GATKDuplicationMetrics extends DuplicationMetrics implements Serializable {
   @NoMergingKeepsValue //Currently picard requires all fields to be annotated.
    private static final long serialVersionUID = 1L;

    public GATKDuplicationMetrics copy() {
        final GATKDuplicationMetrics copy = new GATKDuplicationMetrics();
        copy.LIBRARY = this.LIBRARY;
        copy.UNPAIRED_READS_EXAMINED = this.UNPAIRED_READS_EXAMINED;
        copy.READ_PAIRS_EXAMINED = this.READ_PAIRS_EXAMINED;
        copy.SECONDARY_OR_SUPPLEMENTARY_RDS = this.SECONDARY_OR_SUPPLEMENTARY_RDS;
        copy.UNMAPPED_READS = this.UNMAPPED_READS;
        copy.UNPAIRED_READ_DUPLICATES = this.UNPAIRED_READ_DUPLICATES;
        copy.READ_PAIR_DUPLICATES = this.READ_PAIR_DUPLICATES;
        copy.READ_PAIR_OPTICAL_DUPLICATES = this.READ_PAIR_OPTICAL_DUPLICATES;
        copy.PERCENT_DUPLICATION = this.PERCENT_DUPLICATION;
        copy.ESTIMATED_LIBRARY_SIZE = this.ESTIMATED_LIBRARY_SIZE;

        return copy;
    }

    /**
     * Update metrics given a record or GATKRead
     */
    public void updateMetrics(final SAMRecord rec) {
        update(rec.getReadUnmappedFlag(), rec.isSecondaryOrSupplementary(), ReadUtils.readHasMappedMate(rec), rec.getDuplicateReadFlag() );
    }

    public void updateMetrics(final GATKRead read) {
        update(read.isUnmapped(), read.isSecondaryAlignment() || read.isSupplementaryAlignment(), ReadUtils.readHasMappedMate(read), read.isDuplicate() );
    }

    private void update(final boolean readUnmappedFlag, final boolean secondaryOrSupplementary, final boolean mappedMate, final boolean isDuplicate) {
        if (readUnmappedFlag) {
            ++this.UNMAPPED_READS;
        } else if (secondaryOrSupplementary) {
            ++this.SECONDARY_OR_SUPPLEMENTARY_RDS;
        } else if (!mappedMate) {
            ++this.UNPAIRED_READS_EXAMINED;
        } else {
            ++this.READ_PAIRS_EXAMINED; // will need to be divided by 2 at the end
        }
        // Update the duplication metrics
        if (isDuplicate && !secondaryOrSupplementary) {
            if (!mappedMate) {
                ++this.UNPAIRED_READ_DUPLICATES;
            } else {
                ++this.READ_PAIR_DUPLICATES;// will need to be divided by 2 at the end
            }
        }
    }

}
