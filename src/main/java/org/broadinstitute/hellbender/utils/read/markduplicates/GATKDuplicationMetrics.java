package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.util.Histogram;
import org.broadinstitute.hellbender.utils.Utils;
import picard.sam.DuplicationMetrics;

import java.io.Serializable;

/**
 * Metrics that are calculated during the process of marking duplicates
 * within a stream of SAMRecords.
 */
@SuppressWarnings("serial")
public final class GATKDuplicationMetrics extends DuplicationMetrics implements Serializable {
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

}
