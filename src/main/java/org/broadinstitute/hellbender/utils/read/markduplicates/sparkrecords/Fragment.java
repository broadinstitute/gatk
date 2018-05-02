package org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadEnds;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;

/**
 * Class representing a single read fragment at a particular start location without a mapped mate.
 *
 * This class holds onto as little information as possible in an attempt to prevent excessive serialization of
 * during the processing step of MarkDuplicatesSpark
 */
public class Fragment extends PairedEnds {
    protected transient int key;

    private final int firstStartPosition;
    private final int firstUnclippedStartPosition;
    private final short firstRefIndex;
    private final boolean R1R;

    protected final int score;

    public Fragment(final GATKRead first, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy) {
        super(partitionIndex, first.getName());

        this.firstUnclippedStartPosition = ReadUtils.getStrandedUnclippedStart(first);
        this.firstStartPosition = first.getAssignedStart();
        this.firstRefIndex = (short)ReadUtils.getReferenceIndex(first, header);
        this.score = scoringStrategy.score(first);
        this.R1R = first.isReverseStrand();
        this.key = ReadsKey.hashKeyForFragment(firstUnclippedStartPosition,
                isR1R(),
                firstRefIndex,
                ReadUtils.getLibrary(first, header));
    }

    @Override
    public Type getType() {
      return Type.FRAGMENT;
    }
    @Override
    // NOTE: This is transient and thus may not exist if the object gets serialized
    public int key() {
        return key;
    }
    @Override
    public int getScore() {
      return score;
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
        return name + " " + firstUnclippedStartPosition;
    }
}
