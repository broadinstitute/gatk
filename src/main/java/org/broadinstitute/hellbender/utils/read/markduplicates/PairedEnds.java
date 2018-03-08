package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * Struct-like class to store information about the paired reads for mark duplicates.
 */
public class PairedEnds implements OpticalDuplicateFinder.PhysicalLocation {
  private int firstStartPosition;
  private transient GATKRead first;
  private transient GATKRead second;

  // Information used to detect optical dupes
  private transient short readGroup = -1;
  private transient short tile = -1;
  private transient short x = -1;
  private transient short y = -1;
  private transient short libraryId = -1;
  private final boolean fragment;
  private Integer score;

  private final String name;

  private int firstUnclippedStartPosition = -1;
  private int firstRefIndex = -1;

  private int secondUnclippedStartPosition = -1;
  private int secondRefIndex = -1;

  private boolean R1R = false;
  private boolean R2R = false;

  private final int partitionIndex;


  private PairedEnds(final GATKRead first, final boolean fragment, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy) {
    name = first.getName();
    this.first = first;
    firstUnclippedStartPosition = ReadUtils.getStrandedUnclippedStart(first);
    firstStartPosition = first.getAssignedStart();
    firstRefIndex = ReadUtils.getReferenceIndex(first, header);
    this.partitionIndex = partitionIndex;
    this.fragment = fragment;
    this.score = scoringStrategy.score(first);
    this.R1R = first.isReverseStrand();
  }



  /**
   * special constructor for creating fragments
   * this only includes the necessary data to locate the read, the rest is unecessary because it will appear in the paired bucket
   *
   */
  private PairedEnds(GATKRead read, int partitionIndex) {
    this.firstUnclippedStartPosition = ReadUtils.getStrandedUnclippedStart(read);
    this.first = read;
    this.partitionIndex = partitionIndex;
    this.fragment = true;
    this.firstRefIndex = -1;
    this.name = null;
  }

  public static PairedEnds newFragment(final GATKRead first, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy) {
    return new PairedEnds(first, true, header, partitionIndex, scoringStrategy);
  }

  //todo: why???
  public static PairedEnds newPairWithMissingSecond(final GATKRead first, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy) {
    return new PairedEnds(first, false, header, partitionIndex, scoringStrategy);
  }

  // An optimization for passing around empty read information
  public static PairedEnds empty(GATKRead read, int partitionIndex) {
    return new PairedEnds(read, partitionIndex);
  }

  public static PairedEnds newPair(GATKRead first, GATKRead second, SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy) {
    Utils.nonNull(first);
    Utils.nonNull(second);
    final PairedEnds incomplete = new PairedEnds(first, false, header, partitionIndex, scoringStrategy);

    incomplete.score = incomplete.score + scoringStrategy.score(second);

    if (incomplete.firstUnclippedStartPosition > ReadUtils.getStrandedUnclippedStart(second)) {
      incomplete.firstStartPosition = second.getAssignedStart();
      incomplete.second = incomplete.first;
      incomplete.secondRefIndex = incomplete.firstRefIndex;
      incomplete.secondUnclippedStartPosition = incomplete.firstUnclippedStartPosition;
      incomplete.first = second;
      incomplete.firstRefIndex = ReadUtils.getReferenceIndex(second, header);
      incomplete.firstUnclippedStartPosition = ReadUtils.getStrandedUnclippedStart(second); //TODO don't compute twice
    } else {
      incomplete.second = second;
      incomplete.secondRefIndex = ReadUtils.getReferenceIndex(second, header);
      incomplete.secondUnclippedStartPosition = ReadUtils.getStrandedUnclippedStart(second); //TODO don't compute twice
    }

    // Calculating necessary optical duplicate information

    final GATKRead read1;
    final GATKRead read2;
    if (incomplete.first.isFirstOfPair()){
      read1 = incomplete.first;
      read2 = second;
    } else {
      read1 = second;
      read2 = incomplete.first;
    }

    incomplete.R1R = read1.isReverseStrand();
    incomplete.R2R = read2.isReverseStrand();

    return incomplete;
  }

  public int key(final SAMFileHeader header) {
    firstUnclippedStartPosition = ReadUtils.getStrandedUnclippedStart(first);
    return ReadsKey.hashKeyForPair(header, first, second);
  }

  public int keyForFragment(final SAMFileHeader header) {
    firstUnclippedStartPosition = ReadUtils.getStrandedUnclippedStart(first); //todo
    return ReadsKey.hashKeyForFragment(header, first);
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
    return firstRefIndex == -1;
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
}
