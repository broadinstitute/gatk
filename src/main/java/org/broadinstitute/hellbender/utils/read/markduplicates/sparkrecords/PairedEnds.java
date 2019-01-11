package org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords;

/**
 * Struct-like class to store information about the paired reads for mark duplicates.
 */
public abstract class PairedEnds extends MarkDuplicatesSparkRecord {

    PairedEnds(int partitionIndex, String name) {
        super(partitionIndex, name);
    }

    public abstract short getScore();
    public abstract boolean isRead1ReverseStrand();

    /**
     * This returns a byte summary spanning from 0-5 representing all combinations of single read or two read
     * forward/reverse strand for the first and second read in the Pair. Note, for PCR Duplicates orientation
     * that the 'first' read corresponds to the read that appears first by coordinate sorting order, which is
     * distinct from optical sorting which considers the 'first' read to be the first in pair.
     */
    public abstract byte getOrientationForPCRDuplicates();
}
