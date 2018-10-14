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
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;
import picard.sam.markduplicates.util.ReadEnds;

import java.util.Map;

/**
 * Class representing a pair of reads together with accompanying optical duplicate marking information.
 *
 * This class holds onto as little information as possible in an attempt to prevent excessive serialization of
 * during the processing step of MarkDuplicatesSpark
 */
@DefaultSerializer(Pair.Serializer.class)
public final class Pair extends TransientFieldPhysicalLocation {
    protected transient ReadsKey key;

    private final boolean isRead1ReverseStrand;

    private final boolean isRead2ReverseStrand;
    private final short score;
    private final boolean wasFlipped;

    public Pair(final GATKRead read1, final GATKRead read2, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy, Map<String, Byte> headerLibraryMap) {
        super(partitionIndex, read1.getName());

        final String name1 = read1.getName();
        final String name2 = read2.getName();
        Utils.validate(name1.equals(name2), () -> "Paired reads have different names\n" + name1 + "\n" + name2);

        this.score = (short)(scoringStrategy.score(read1) + scoringStrategy.score(read2));

        GATKRead first = read1;
        GATKRead second;

        final int read1UnclippedStart = ReadUtils.getStrandedUnclippedStart(read1);
        final int read2UnclippedStart = ReadUtils.getStrandedUnclippedStart(read2);

        int read1ReferenceIndex = ReadUtils.getReferenceIndex(read1,header);
        int read2ReferenceIndex = ReadUtils.getReferenceIndex(read2,header);

        if( read1ReferenceIndex != read2ReferenceIndex ? read1ReferenceIndex < read2ReferenceIndex : read1UnclippedStart <= read2UnclippedStart ){
            first = read1;
            second = read2;
        } else {
            first = read2;
            second = read1;
        }

        // if the two read ends are in the same position, pointing in opposite directions,
        // the orientation is undefined and the procedure above
        // will depend on the order of the reads in the file.
        // To avoid this, and match picard's behavior in this case, we ensure the orientation will be FR:
        if (read1ReferenceIndex == read2ReferenceIndex &&
                read1UnclippedStart == read2UnclippedStart &&
                first.isReverseStrand() && !second.isReverseStrand()) {
            // setting to FR for consistencies sake. (which involves flipping) if both reads had the same unclipped start
            GATKRead tmp = first;
            first = second;
            second = tmp;
        }

        this.key = ReadsKey.getKeyForPair(header, first, second, headerLibraryMap);

        isRead1ReverseStrand = first.isReverseStrand();
        isRead2ReverseStrand = second.isReverseStrand();

        // Keep track of the orientation of read1 and read2 as it is important for optical duplicate marking
        wasFlipped = second.isFirstOfPair();

        this.key = ReadsKey.getKeyForPair(header, first, second, headerLibraryMap);
    }

    // Constructor for serialization purposes
    private Pair(Kryo kryo, Input input){
        super(input.readInt(true), input.readString());

        // Information used to detect optical dupes
        tile = -1;
        x = -1;
        y = -1;
        libraryId = -1;

        score = input.readShort();

        isRead1ReverseStrand = input.readBoolean();

        isRead2ReverseStrand = input.readBoolean();

        readGroupIndex = input.readShort();
        wasFlipped = input.readBoolean();
    }

    protected void serialize(Kryo kryo, Output output) {
        output.writeInt(partitionIndex, true);
        output.writeAscii(name);

        output.writeShort(score);

        output.writeBoolean(isRead1ReverseStrand);

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
    public ReadsKey key() {
        return key;
    }
    @Override
    public short getScore() {
        return score;
    }

    @Override
    public boolean isRead1ReverseStrand() {
        return isRead1ReverseStrand;
    }
    @Override
    public String toString() {
        return name + " score:" + score;
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
