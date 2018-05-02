package org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadEnds;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;

/**
 * Class representing a pair of reads together with accompanying optical duplicate marking information.
 *
 * This class holds onto as little information as possible in an attempt to prevent excessive serialization of
 * during the processing step of MarkDuplicatesSpark
 */
@DefaultSerializer(Pair.Serializer.class)
public final class Pair extends PairedEnds implements picard.sam.util.PhysicalLocation {
    protected transient int key;

    private final int firstStartPosition;
    private final int firstUnclippedStartPosition;
    private final short firstRefIndex;
    private final boolean R1R;

    private final int secondUnclippedStartPosition;
    private final short secondRefIndex;
    private final boolean R2R;
    private final int score;

    // Information used to detect optical dupes
    private short readGroupIndex = -1;
    private transient short tile = -1;
    private transient short x = -1;
    private transient short y = -1;
    private transient short libraryId = -1;

    public Pair(final GATKRead read1, final GATKRead read2, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy) {
        super(partitionIndex, read1.getName());

        final String name1 = read1.getName();
        final String name2 = read2.getName();
        Utils.validate(name1.equals(name2), () -> "Paired reads have different names\n" + name1 + "\n" + name2);

        this.score = scoringStrategy.score(read1) + scoringStrategy.score(read2);

        GATKRead first = read1;
        GATKRead second;

        final int read1UnclippedStart = ReadUtils.getStrandedUnclippedStart(read1);
        final int read2UnclippedStart = ReadUtils.getStrandedUnclippedStart(read2);

        if( read1UnclippedStart < read2UnclippedStart ){
            first = read1;
            second = read2;
            firstUnclippedStartPosition = read1UnclippedStart;
            secondUnclippedStartPosition = read2UnclippedStart;
        } else {
            first = read2;
            second = read1;
            firstUnclippedStartPosition = read2UnclippedStart;
            secondUnclippedStartPosition = read1UnclippedStart;
        }

        firstStartPosition = first.getAssignedStart();
        firstRefIndex = (short)ReadUtils.getReferenceIndex(first, header);
        secondRefIndex = (short)ReadUtils.getReferenceIndex(second, header);

        final GATKRead firstOfPair;
        final GATKRead secondOfPair;
        if (read1.isFirstOfPair()){
            firstOfPair = read1;
            secondOfPair = read2;
        } else {
            firstOfPair = read2;
            secondOfPair = read1;
        }

        R1R = firstOfPair.isReverseStrand();
        R2R = secondOfPair.isReverseStrand();

        this.key = ReadsKey.hashKeyForPair(header, first, second);
    }

    // Constructor for serialization purposes
    private Pair(Kryo kryo, Input input){
        super(input.readInt(true), input.readString());

        // Information used to detect optical dupes
        readGroupIndex = -1;
        tile = -1;
        x = -1;
        y = -1;
        libraryId = -1;

        score = input.readInt();

        firstStartPosition = input.readInt();
        firstUnclippedStartPosition = input.readInt();
        firstRefIndex = input.readShort();
        R1R = input.readBoolean();

        secondUnclippedStartPosition = input.readInt();
        secondRefIndex = input.readShort();
        R2R = input.readBoolean();

        readGroupIndex = input.readShort();
    }

    protected void serialize(Kryo kryo, Output output) {
        output.writeInt(partitionIndex, true);
        output.writeAscii(name);

        output.writeInt(score);

        output.writeInt(firstStartPosition);
        output.writeInt(firstUnclippedStartPosition);
        output.writeShort(firstRefIndex);
        output.writeBoolean(R1R);

        output.writeInt(secondUnclippedStartPosition);
        output.writeShort(secondRefIndex);
        output.writeBoolean(R2R);

        output.writeShort(readGroupIndex);
    }

    @Override
    public Type getType() {
        return Type.PAIR;
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
    public int getFirstRefIndex() {
        return firstRefIndex;
    }
    @Override
    public String toString() {
        return name + " " + firstUnclippedStartPosition + " " + (secondUnclippedStartPosition == -1 ? "" : secondUnclippedStartPosition);
    }

    /**
     * Returns the pair orientation suitable for optical duplicates,
     * which always goes by the first then the second end for the strands.
     * This is based on code in MarkDuplicatesGATK and ReadEnds.getOrientationByte.
     * Returns one of {@link ReadEnds#RR}, {@link ReadEnds#RF}, {@link ReadEnds#FR}, {@link ReadEnds#FF}
     */
    public byte getOrientationForOpticalDuplicates() {
        if (R1R && R2R) {
            return ReadEnds.RR;
        }
        if (R1R) {
            return ReadEnds.RF; //at this point we know for sure R2R is false
        }
        if (R2R) {
            return ReadEnds.FR; //at this point we know for sure R1R is false
        }
        return ReadEnds.FF;  //at this point we know for sure R1R is false and R2R is false
    }

    // Methods for OpticalDuplicateFinder.PhysicalLocation
    @Override
    public short getReadGroup() { return this.readGroupIndex; }
    @Override
    public void setReadGroup(final short readGroup) { this.readGroupIndex = readGroup; }
    @Override
    public short getTile() { return this.tile; }
    @Override
    public void setTile(final short tile) { this.tile = tile; }
    @Override
    public short getX() { return this.x; }
    @Override
    public void setX(final short x) { this.x = x; }
    @Override
    public short getY() { return this.y; }
    @Override
    public void setY(final short y) { this.y = y; }
    @Override
    public short getLibraryId() { return this.libraryId; }
    @Override
    public void setLibraryId(final short libraryId) { this.libraryId = libraryId; }

    /**
     * Serializers for each subclass of PairedEnds which rely on implementations of serializations within each class itself
     */
    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords.Pair> {
        @Override
        public void write(final Kryo kryo, final Output output, final org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords.Pair pair ) {
            pair.serialize(kryo, output);
        }
        @Override
        public org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords.Pair read(final Kryo kryo, final Input input, final Class<org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords.Pair> klass ) {
            return new org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords.Pair(kryo, input);
        }
    }
}
