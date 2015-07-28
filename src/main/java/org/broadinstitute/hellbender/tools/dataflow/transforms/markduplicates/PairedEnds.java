package org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * Struct-like class to store information about the paired reads for mark duplicates.
 */
public final class PairedEnds {
    private GATKRead first, second;

    PairedEnds(final GATKRead first) {
        this.first = first;
    }

    public static PairedEnds of(final GATKRead first) {
        return new PairedEnds(first);
    }

    public PairedEnds and(final GATKRead second) {
        if (second != null &&
                ReadUtils.getStrandedUnclippedStart(first) > ReadUtils.getStrandedUnclippedStart(second)) {

            this.second = this.first;
            this.first = second;
        } else {
            this.second = second;
        }
        return this;
    }

    public String key(final SAMFileHeader header) {
        return ReadsKey.keyForPairedEnds(header, first, second);
    }

    public GATKRead first() {
        return first;
    }

    public GATKRead second() {
        return second;
    }

    public int score() {
        return MarkDuplicatesUtils.score(first) + MarkDuplicatesUtils.score(second);
    }
}
