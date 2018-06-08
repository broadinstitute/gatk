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
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadEnds;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;
import picard.sam.util.PhysicalLocation;

/**
 * Class representing a pair of reads together with accompanying optical duplicate marking information.
 *
 * This class holds onto as little information as possible in an attempt to prevent excessive serialization of
 * during the processing step of MarkDuplicatesSpark
 */
@DefaultSerializer(Pair.Serializer.class)
public final class Pair extends PairedEnds implements PhysicalLocation {
    protected transient String key;

    private final int firstStartPosition;
    private final int firstUnclippedStartPosition;
    private final short firstRefIndex;
    private final boolean isRead1ReverseStrand;

    private final int secondUnclippedStartPosition;
    private final short secondRefIndex;
    private final boolean isRead2ReverseStrand;
    private final int score;
    private final boolean wasFlipped;

    // Information used to detect optical dupes
    private short readGroupIndex = -1;

    private transient short tile = -1;
    private transient int x = -1;
    private transient int y = -1;
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
        // Keep track of the orientation of read1 and read2 as it is important for optical duplicate marking
        wasFlipped = second.isFirstOfPair();

        firstStartPosition = first.getAssignedStart();
        firstRefIndex = (short)ReadUtils.getReferenceIndex(first, header);
        secondRefIndex = (short)ReadUtils.getReferenceIndex(second, header);

        isRead1ReverseStrand = first.isReverseStrand();
        isRead2ReverseStrand = second.isReverseStrand();

        this.key = ReadsKey.hashKeyForPair(header, first, second);
    }

    // Constructor for serialization purposes
    private Pair(Kryo kryo, Input input){
        super(input.readInt(true), input.readString());

        // Information used to detect optical dupes
        tile = -1;
        x = -1;
        y = -1;
        libraryId = -1;

        score = input.readInt();

        firstStartPosition = input.readInt();
        firstUnclippedStartPosition = input.readInt();
        firstRefIndex = input.readShort();
        isRead1ReverseStrand = input.readBoolean();

        secondUnclippedStartPosition = input.readInt();
        secondRefIndex = input.readShort();
        isRead2ReverseStrand = input.readBoolean();

        readGroupIndex = input.readShort();
        wasFlipped = input.readBoolean();
    }

    protected void serialize(Kryo kryo, Output output) {
        output.writeInt(partitionIndex, true);
        output.writeAscii(name);

        output.writeInt(score);

        output.writeInt(firstStartPosition);
        output.writeInt(firstUnclippedStartPosition);
        output.writeShort(firstRefIndex);
        output.writeBoolean(isRead1ReverseStrand);

        output.writeInt(secondUnclippedStartPosition);
        output.writeShort(secondRefIndex);
        output.writeBoolean(isRead2ReverseStrand);

        output.writeShort(readGroupIndex);
        output.writeBoolean(wasFlipped);
    }

    @Override
    public Type getType() {
        return Type.PAIR;
    }
    @Override
    // NOTE: This is transient and thus may not exist if the object gets serialized
    public String key() {
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
    public boolean isRead1ReverseStrand() {
        return isRead1ReverseStrand;
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
        return getOrientation(true);
    }

    @Override
    public byte getOrientationForPCRDuplicates() {
        return getOrientation(false);
    }

    private byte getOrientation(final boolean optical) {
        if (isRead1ReverseStrand && isRead2ReverseStrand) {
            return ReadEnds.RR;
        }
        if (isRead1ReverseStrand) {
            return (optical && wasFlipped)? ReadEnds.FR : ReadEnds.RF; //at this point we know for sure isRead2ReverseStrand is false
        }
        if (isRead2ReverseStrand) {
            return (optical && wasFlipped)? ReadEnds.RF :ReadEnds.FR; //at this point we know for sure isRead1ReverseStrand is false
        }
        return ReadEnds.FF;  //at this point we know for sure isRead1ReverseStrand is false and isRead2ReverseStrand is false
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
    public int getX() { return this.x; }
    @Override
    public void setX(final int x) { this.x = x; }
    @Override
    public int getY() { return this.y; }
    @Override
    public void setY(final int y) { this.y = y; }
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
