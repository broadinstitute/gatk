package org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.LibraryIdGenerator;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;
import picard.sam.markduplicates.util.ReadEnds;

import java.util.Map;

/**
 * Dummy class representing a mated read fragment at a particular start position to be used for accounting
 * when deciding whether to duplicate unmatched fragments.
 *
 * This class holds onto as little information as possible in an attempt to prevent excessive serialization of
 * during the processing step of MarkDuplicatesSpark
 */
public final class EmptyFragment extends PairedEnds {
    protected transient ReadsKey key;

    private final boolean R1R;

    /**
     * special constructor for creating empty fragments
     * this only includes the necessary data to locate the read, the rest is unnecessary because it will appear in the paired bucket
     *
     */
    public EmptyFragment(GATKRead read, SAMFileHeader header, Map<String, Byte> headerLibraryMap) {
        super(0, null);
        this.R1R = read.isReverseStrand();
        this.key = ReadsKey.getKeyForFragment(ReadUtils.getStrandedUnclippedStart(read),
                isRead1ReverseStrand(),
                ReadUtils.getReferenceIndex(read, header),
                headerLibraryMap.get(MarkDuplicatesSparkUtils.getLibraryForRead(read, header, LibraryIdGenerator.UNKNOWN_LIBRARY)));
    }

    @Override
    public Type getType() {
        return Type.EMPTY_FRAGMENT;
    }
    @Override
    public short getScore() {
        return 0;
    }
    @Override
    // NOTE: This is transient and thus may not exist if the object gets serialized
    public ReadsKey key() {
        return key;
    }
    @Override
    public boolean isRead1ReverseStrand() {
        return R1R;
    }
    @Override
    public byte getOrientationForPCRDuplicates() {
        return (R1R)? ReadEnds.R : ReadEnds.F;
    }
    @Override
    public String toString() {
        return "EmptyFragment ";
    }
}
