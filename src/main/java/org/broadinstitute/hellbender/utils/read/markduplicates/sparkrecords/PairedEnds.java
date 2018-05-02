package org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords;

/**
 * Struct-like class to store information about the paired reads for mark duplicates.
 */
public abstract class PairedEnds extends MarkDuplicatesSparkRecord {

    PairedEnds(int partitionIndex, String name) {
        super(partitionIndex, name);
    }

    public abstract int getFirstStartPosition();
    public abstract int getUnclippedStartPosition();
    public abstract int getFirstRefIndex();
    public abstract byte getPCROrientation();
    public abstract int getScore();
    public abstract boolean isR1R();
}
