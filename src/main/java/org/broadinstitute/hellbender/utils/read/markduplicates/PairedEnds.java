package org.broadinstitute.hellbender.utils.read.markduplicates;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * Struct-like class to store information about the paired reads for mark duplicates.
 */
public abstract class PairedEnds extends MarkDuplicatesSparkData {

    PairedEnds(int partitionIndex, String name) {
        super(partitionIndex, name);
    }

    public static PairedEnds newFragment(final GATKRead first, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy) {
        return new Fragment(first, header, partitionIndex, scoringStrategy);
    }

    // An optimization for passing around empty read information
    public static PairedEnds placeHolder(GATKRead read, SAMFileHeader header, int partitionIndex) {
        return new EmptyFragment(read, header, partitionIndex);
    }

    public static PairedEnds newPair(GATKRead first, GATKRead second, SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy) {
        return new Pair(first, second, header, partitionIndex, scoringStrategy);
    }

    public abstract int getFirstStartPosition();
    public abstract boolean isR1R();

    /**
     * Class representing a pair of reads together with accompanying optical duplicate marking information.
     *
     * This class holds onto as little information as possible in an attempt to prevent excessive serialization of
     * during the processing step of MarkDuplicatesSpark
     */
    @DefaultSerializer(PairedEnds.PairedEndsPairSerializer.class)
    public static final class Pair extends PairedEnds implements OpticalDuplicateFinder.PhysicalLocation {
        protected transient GATKRead first;
        protected transient GATKRead second;

        private final int firstStartPosition;
        private final int firstUnclippedStartPosition;
        private final short firstRefIndex;
        private final boolean R1R;

        private final int secondUnclippedStartPosition;
        private final short secondRefIndex;
        private final boolean R2R;
        private final int score;

        // Information used to detect optical dupes
        private transient short readGroup = -1;
        private transient short tile = -1;
        private transient short x = -1;
        private transient short y = -1;
        private transient short libraryId = -1;

        // Constructor for serialization purposes
        private Pair(final GATKRead read1, final GATKRead read2, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy) {
            super(partitionIndex, read1.getName());

            final String name1 = read1.getName();
            final String name2 = read2.getName();
            Utils.validate(name1.equals(name2), () -> "Paired reads have different names\n" + name1 + "\n" + name2);

            this.score = scoringStrategy.score(read1) + scoringStrategy.score(read2);

            this.first = read1;
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
        }

        // Constructor for serialization purposes
        private Pair(Kryo kryo, Input input){
            super(input.readInt(true), input.readString());
            first = null;
            second = null;

            // Information used to detect optical dupes
            readGroup = -1;
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
        }

        @Override
        public Type getType() {
            return Type.PAIR;
        }
        @Override
        public int key(SAMFileHeader header) {
            return ReadsKey.hashKeyForPair(header, first, second);
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
        public short getReadGroup() { return this.readGroup; }
        @Override
        public void setReadGroup(final short readGroup) { this.readGroup = readGroup; }
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
    }


    /**
     * Class representing a single read fragment at a particular start location without a mapped mate.
     *
     * This class holds onto as little information as possible in an attempt to prevent excessive serialization of
     * during the processing step of MarkDuplicatesSpark
     */
    @DefaultSerializer(PairedEnds.PairedEndsFragmentSerializer.class)
    public static class Fragment extends PairedEnds{
        protected transient GATKRead first;

        private final int firstStartPosition;
        private final int firstUnclippedStartPosition;
        private final short firstRefIndex;
        private final boolean R1R;

        protected final int score;

        private Fragment(final GATKRead first, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy) {
            super(partitionIndex, first.getName());

            this.first = first;
            this.firstUnclippedStartPosition = ReadUtils.getStrandedUnclippedStart(first);
            this.firstStartPosition = first.getAssignedStart();
            this.firstRefIndex = (short)ReadUtils.getReferenceIndex(first, header);
            this.score = scoringStrategy.score(first);
            this.R1R = first.isReverseStrand();
        }

        // Constructor for serialization purposes
        private Fragment(Kryo kryo, Input input){
            super(input.readInt(true), input.readString());
            this.first = null;

            this.score = input.readInt();

            this.firstStartPosition = input.readInt();
            this.firstUnclippedStartPosition = input.readInt();
            this.firstRefIndex = input.readShort();
            this.R1R = input.readBoolean();
        }
        @Override
        protected void serialize(Kryo kryo, Output output) {
            output.writeInt(partitionIndex, true);
            output.writeAscii(name);

            output.writeInt(score);

            output.writeInt(firstStartPosition);
            output.writeInt(firstUnclippedStartPosition);
            output.writeShort(firstRefIndex);
            output.writeBoolean(R1R);
        }

        @Override
        public Type getType() {
          return Type.FRAGMENT;
        }
        @Override
        public int key(SAMFileHeader header) {
            return ReadsKey.hashKeyForFragment(firstUnclippedStartPosition,
                    isR1R(),
                    firstRefIndex,
                    ReadUtils.getLibrary(first, header));
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
            return name + " " + firstUnclippedStartPosition;
        }
    }

    /**
     * Dummy class representing a mated read fragment at a particular start position to be used for accounting
     * when deciding whether to duplicate unmatched fragments.
     *
     * This class holds onto as little information as possible in an attempt to prevent excessive serialization of
     * during the processing step of MarkDuplicatesSpark
     */
    @DefaultSerializer(PairedEnds.PairedEndsEmptyFragmentSerializer.class)
    public static final class EmptyFragment extends PairedEnds{
        protected transient GATKRead first;

        private final int firstUnclippedStartPosition;
        private final int firstStartPosition;
        private final short firstRefIndex;
        private final boolean R1R;

        /**
         * special constructor for creating empty fragments
         * this only includes the necessary data to locate the read, the rest is unnecessary because it will appear in the paired bucket
         *
         */
        private EmptyFragment(GATKRead read, SAMFileHeader header, int partitionIndex) {
            super(0, null);

            this.firstUnclippedStartPosition = ReadUtils.getStrandedUnclippedStart(read);
            this.first = read;
            this.firstRefIndex = (short)ReadUtils.getReferenceIndex(read, header);
            this.R1R = read.isReverseStrand();
            firstStartPosition = 0;
        }

        // Constructor for serialization purposes
        private EmptyFragment(Kryo kryo, Input input){
            super(input.readInt(true), null);
            first = null;

            firstStartPosition = input.readInt();
            firstUnclippedStartPosition = input.readInt();
            firstRefIndex = input.readShort();
            R1R = input.readBoolean();
        }
        @Override
        protected void serialize(Kryo kryo, Output output) {
            output.writeInt(partitionIndex, true);

            output.writeInt(firstStartPosition);
            output.writeInt(firstUnclippedStartPosition);
            output.writeShort(firstRefIndex);
            output.writeBoolean(R1R);
        }

        @Override
        public Type getType() {
            return Type.FRAGMENT;
        }
        @Override
        public int getScore() {
            return 0;
        }
        @Override
        public int key(SAMFileHeader header) {
            return ReadsKey.hashKeyForFragment(firstUnclippedStartPosition,
                    isR1R(),
                    firstRefIndex,
                    ReadUtils.getLibrary(first, header));
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
            return "EmptyFragment " + firstUnclippedStartPosition;
        }
    }

    //===============================================================================================================================
    // Serializers
    //===============================================================================================================================
    /**
     * Serializers for each subclass of PairedEnds which rely on implementations of serializations within each class itself
     */
    public static final class PairedEndsPairSerializer extends com.esotericsoftware.kryo.Serializer<PairedEnds.Pair> {
        @Override
        public void write(final Kryo kryo, final Output output, final PairedEnds.Pair pair ) {
            pair.serialize(kryo, output);
        }
        @Override
        public PairedEnds.Pair read(final Kryo kryo, final Input input, final Class<PairedEnds.Pair> klass ) {
            return new PairedEnds.Pair(kryo, input);
        }
    }
    public static final class PairedEndsFragmentSerializer extends com.esotericsoftware.kryo.Serializer<PairedEnds.Fragment> {
        @Override
        public void write(final Kryo kryo, final Output output, final PairedEnds.Fragment pair ) {
            pair.serialize(kryo, output);
        }
        @Override
        public PairedEnds.Fragment read(final Kryo kryo, final Input input, final Class<PairedEnds.Fragment> klass ) {
            return new PairedEnds.Fragment(kryo, input);
        }
    }
    public static final class PairedEndsEmptyFragmentSerializer extends com.esotericsoftware.kryo.Serializer<PairedEnds.EmptyFragment> {
        @Override
        public void write(final Kryo kryo, final Output output, final PairedEnds.EmptyFragment pair ) {
            pair.serialize(kryo, output);
        }
        @Override
        public PairedEnds.EmptyFragment read(final Kryo kryo, final Input input, final Class<PairedEnds.EmptyFragment> klass ) {
            return new PairedEnds.EmptyFragment(kryo, input);
        }
    }
}
