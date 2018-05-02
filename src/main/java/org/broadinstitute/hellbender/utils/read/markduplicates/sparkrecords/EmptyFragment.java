package org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadEnds;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;

/**
 * Dummy class representing a mated read fragment at a particular start position to be used for accounting
 * when deciding whether to duplicate unmatched fragments.
 *
 * This class holds onto as little information as possible in an attempt to prevent excessive serialization of
 * during the processing step of MarkDuplicatesSpark
 */
public final class EmptyFragment extends PairedEnds {
    protected transient int key;

    private final int firstUnclippedStartPosition;
    private final int firstStartPosition;
    private final short firstRefIndex;
    private final boolean R1R;

    /**
     * special constructor for creating empty fragments
     * this only includes the necessary data to locate the read, the rest is unnecessary because it will appear in the paired bucket
     *
     */
    public EmptyFragment(GATKRead read, SAMFileHeader header) {
        super(0, null);

        this.firstUnclippedStartPosition = ReadUtils.getStrandedUnclippedStart(read);
        this.firstRefIndex = (short)ReadUtils.getReferenceIndex(read, header);
        this.R1R = read.isReverseStrand();
        firstStartPosition = 0;
        this.key = ReadsKey.hashKeyForFragment(firstUnclippedStartPosition,
                isR1R(),
                firstRefIndex,
                ReadUtils.getLibrary(read, header));
    }

    @Override
    public Type getType() {
        return Type.EMPTY_FRAGMENT;
    }
    @Override
    public int getScore() {
        return 0;
    }
    @Override
    // NOTE: This is transient and thus may not exist if the object gets serialized
    public int key() {
        return key;
    }
    @Override
    public int getUnclippedStartPosition() {
        return firstUnclippedStartPosition;
    }
    @Override
    public int getFirstStartPosition() {
        return firstStartPosition;
    }
    @Override
    public boolean isR1R() {
        return R1R;
    }
    @Override
    public byte getPCROrientation() {
        return (R1R)? ReadEnds.R : ReadEnds.F;
    }
    @Override
    public int getFirstRefIndex() {
        return firstRefIndex;
    }
    @Override
    public String toString() {
        return "EmptyFragment " + firstUnclippedStartPosition;
    }
}
