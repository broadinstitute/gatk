package org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;

/**
 * Dummy class used for preserving reads that need to be marked as non-duplicate despite not wanting to perform any
 * processing on the reads. (eg. unmapped reads we don't want to process but must be non-duplicate marked)
 */
public final class Passthrough extends MarkDuplicatesSparkRecords {
    private final transient GATKRead read;

    Passthrough(GATKRead read, int partitionIndex) {
        super(partitionIndex, read.getName());

        this.read = read;
    }

    @Override
    public Type getType() {
        return Type.PASSTHROUGH;
    }
    @Override
    public int key(final SAMFileHeader header) {
        return ReadsKey.hashKeyForRead(read);
    }
}
