package org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords;

import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;

/**
 * Dummy class used for preserving reads that need to be marked as non-duplicate despite not wanting to perform any
 * processing on the reads. (eg. unmapped reads we don't want to process but must be non-duplicate marked)
 */
public final class Passthrough extends MarkDuplicatesSparkRecord {
    private final transient ReadsKey key;

    Passthrough(GATKRead read, int partitionIndex) {
        super(partitionIndex, read.getName());

        // use a hash key here instead of a normal key because collisions don't matter here
        this.key = ReadsKey.hashKeyForPassthroughRead(read);
    }

    @Override
    public Type getType() {
        return Type.PASSTHROUGH;
    }
    @Override
    // NOTE: This is transient and thus may not exist if the object gets serialized
    public ReadsKey key() {
        return key;
    }
}
