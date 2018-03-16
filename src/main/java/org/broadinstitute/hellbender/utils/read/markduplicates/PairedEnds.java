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
@DefaultSerializer(PairedEnds.Serializer.class)
public class PairedEnds implements OpticalDuplicateFinder.PhysicalLocation {

  private transient GATKRead first;
  private transient GATKRead second;

  // Information used to detect optical dupes
  private transient short readGroup = -1;
  private transient short tile = -1;
  private transient short x = -1;
  private transient short y = -1;
  private transient short libraryId = -1;

  private final int partitionIndex;
  private final boolean fragment;
  private final int score;
  private final String name;

  private final int firstStartPosition;
  private final int firstUnclippedStartPosition;
  private final short firstRefIndex;
  private final boolean R1R;

  private final int secondUnclippedStartPosition;
  private final short secondRefIndex;
  private final boolean R2R;



  private PairedEnds(final GATKRead first, final boolean fragment, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy) {
    this.name = first.getName();
    this.first = first;
    this.firstUnclippedStartPosition = ReadUtils.getStrandedUnclippedStart(first);
    this.secondUnclippedStartPosition = -1;
    this.firstStartPosition = first.getAssignedStart();
    this.firstRefIndex = (short)ReadUtils.getReferenceIndex(first, header);
    this.secondRefIndex = -1;
    this.partitionIndex = partitionIndex;
    this.fragment = fragment;
    this.score = scoringStrategy.score(first);
    this.R1R = first.isReverseStrand();
    this.R2R = false;
  }


  /**
   * special constructor for creating empty fragments
   * this only includes the necessary data to locate the read, the rest is unnecessary because it will appear in the paired bucket
   *
   */
  private PairedEnds(GATKRead read, SAMFileHeader header, int partitionIndex) {
    this.firstUnclippedStartPosition = ReadUtils.getStrandedUnclippedStart(read);
    this.secondUnclippedStartPosition = -1;
    this.first = read;
    this.partitionIndex = partitionIndex;
    this.fragment = true;
    this.firstRefIndex = (short)ReadUtils.getReferenceIndex(read, header);
    this.secondRefIndex = -1;
    this.name = null;
    this.R1R = read.isReverseStrand();
    this.R2R = false;
    this.score = 0;
    this.firstStartPosition = firstUnclippedStartPosition;
  }

  private PairedEnds(final GATKRead read1, final GATKRead read2, final boolean fragment, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy) {
    final String name1 = read1.getName();
    final String name2 = read2.getName();
    Utils.validate(name1.equals(name2), () -> "Paired reads have different names\n" + name1 + "\n" + name2);
    this.name = name1;

    this.partitionIndex = partitionIndex;
    this.fragment = fragment;
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

  public static PairedEnds newFragment(final GATKRead first, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy) {
    return new PairedEnds(first, true, header, partitionIndex, scoringStrategy);
  }

  //todo: why???
  public static PairedEnds newPairWithMissingSecond(final GATKRead first, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy) {
    return new PairedEnds(first, false, header, partitionIndex, scoringStrategy);
  }

  // An optimization for passing around empty read information
  public static PairedEnds placeHolder(GATKRead read, SAMFileHeader header, int partitionIndex) {
    return new PairedEnds(read, header, partitionIndex);
  }

  public static PairedEnds newPair(GATKRead first, GATKRead second, SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy) {
    return new PairedEnds(first, second, false, header, partitionIndex, scoringStrategy);
  }

  public Type getType(){
      if(isFragment()){
          return Type.FRAGMENT;
      } else if (hasSecondRead()) {
          return Type.PAIR;
      } else {
          return Type.PAIRED_BUT_MISSING_SECOND_READ;
      }
  }

  public int key(final SAMFileHeader header) {
    return ReadsKey.hashKeyForPair(header, first, second);
  }

  public int keyForFragment(final SAMFileHeader header) {
    return ReadsKey.hashKeyForFragment(firstUnclippedStartPosition,
                                       isR1R(),
                                       firstRefIndex,
                                       ReadUtils.getLibrary(first, header));
  }

  public int getScore() {
    return score;
  }

  public int getUnclippedStartPosition() {
    return firstUnclippedStartPosition;
  }

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

  public boolean isFragment() {
    return fragment;
  }

  //TODO note that this only
  public boolean isEmpty() {
    return name == null;
  }

  /**
   * Returns the pair orientation suitable for optical duplicates,
   * which always goes by the first then the second end for the strands.
   * This is based on code in MarkDuplicatesGATK and ReadEnds.getOrientationByte.
   * Returns one of {@link ReadEnds#RR}, {@link ReadEnds#RF}, {@link ReadEnds#FR}, {@link ReadEnds#FF}
   */
  public byte getOrientationForOpticalDuplicates() {
//    final GATKRead read1;
//    final GATKRead read2;
//    if (first.isFirstOfPair()){
//      read1 = first;
//      read2 = second;
//    } else {
//      read1 = second;
//      read2 = first;
//    }
//
//    final boolean R1R = read1.isReverseStrand();
//    final boolean R2R = read2.isReverseStrand();
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

  public int getPartitionIndex(){
    return partitionIndex;
  }

  public boolean hasSecondRead(){
    return secondRefIndex != -1;
  }

  public String getName() {
    return name;
  }

  public int getFirstStartPosition() {
    return firstStartPosition;
  }

  public boolean isR1R() {
    return R1R;
  }

  public int getFirstRefIndex() {
    return firstRefIndex;
  }

  @Override
  public String toString() {
    return name + " " + firstUnclippedStartPosition + " " + (secondUnclippedStartPosition == -1 ? "" : secondUnclippedStartPosition);
  }

  public enum Type {
      FRAGMENT, PAIR, PAIRED_BUT_MISSING_SECOND_READ;

  }

  public static final class Serializer extends com.esotericsoftware.kryo.Serializer<PairedEnds> {
    @Override
    public void write(final Kryo kryo, final Output output, final PairedEnds pair ) {
      pair.serialize(kryo, output);
    }

    @Override
    public PairedEnds read(final Kryo kryo, final Input input, final Class<PairedEnds> klass ) {
      return new PairedEnds(kryo, input);
    }
  }


  private void serialize(Kryo kryo, Output output) {
//    private final int partitionIndex;
//    private final boolean fragment;
//    private final int score;
//    private final String name;
//
//    private final int firstStartPosition;
//    private final int firstUnclippedStartPosition;
//    private final short firstRefIndex;
//    private final boolean R1R;
//
//    private final int secondUnclippedStartPosition;
//    private final short secondRefIndex;
//    private final boolean R2R;
    
    output.writeInt(partitionIndex, true);
    output.writeBoolean(fragment);
    output.writeInt(score);
    output.writeAscii(name);

    output.writeInt(firstStartPosition);
    output.writeInt(firstUnclippedStartPosition);
    output.writeShort(firstRefIndex);
    output.writeBoolean(R1R);

    if(!isFragment()){
      output.writeInt(secondUnclippedStartPosition);
      output.writeShort(secondRefIndex);
      output.writeBoolean(R2R);
    }

  }

  private PairedEnds(Kryo kryo, Input input){

    first = null;
    second = null;

    // Information used to detect optical dupes
    readGroup = -1;
    tile = -1;
    x = -1;
    y = -1;
    libraryId = -1;

    partitionIndex = input.readInt(true);
    fragment = input.readBoolean();
    score = input.readInt();
    name = input.readString();

    firstStartPosition = input.readInt();
    firstUnclippedStartPosition = input.readInt();
    firstRefIndex = input.readShort();
    R1R = input.readBoolean();

    if(!isFragment()){
           secondUnclippedStartPosition = input.readInt();
           secondRefIndex = input.readShort();
           R2R = input.readBoolean();
    } else {
        secondUnclippedStartPosition = -1;
        secondRefIndex = -1;
        R2R = false;
    }
  }

}
