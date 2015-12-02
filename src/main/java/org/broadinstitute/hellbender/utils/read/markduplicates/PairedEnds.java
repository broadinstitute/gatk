package org.broadinstitute.hellbender.utils.read.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

/**
 * Struct-like class to store information about the paired reads for mark duplicates.
 */
public class PairedEnds implements OpticalDuplicateFinder.PhysicalLocation {
  private GATKRead first, second;

  // Information used to detect optical dupes
  public short readGroup = -1;
  public short tile = -1;
  public short x = -1, y = -1;
  public short libraryId = -1;

  PairedEnds(final GATKRead first) {
    this.first = first;
  }

  public static PairedEnds of(final GATKRead first) {
    return new PairedEnds(first);
  }

  public PairedEnds and(final GATKRead second) {
    if (second != null &&
        ReadUtils.getStrandedUnclippedStart(first) > ReadUtils.getStrandedUnclippedStart(second)) {
      this.second = this.first;
      this.first = second;
    } else {
      this.second = second;
    }
    return this;
  }

  public String key(final SAMFileHeader header) {
    return ReadsKey.keyForPairedEnds(header, first, second);
  }

  public String keyForFragment(final SAMFileHeader header) {
    return ReadsKey.keyForFragment(header, first);
  }

  public GATKRead first() {
    return first;
  }

  public GATKRead second() {
    return second;
  }

  public int score(final MarkDuplicatesScoringStrategy scoringStrategy) {
    return scoringStrategy.score(first) + scoringStrategy.score(second);
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
  /**
   * returns a deep(ish) copy of the GATK reads in the PairedEnds.
   * TODO: This is only deep for the Google Model read, GATKRead copy() isn't deep for
   * TODO: for the SAMRecord backed read.
   * @return a new deep copy
   */
  public PairedEnds copy() {
    if (second == null) {
      return new PairedEnds(first.copy());
    }
    return new PairedEnds(first.copy()).and(second.copy());
  }

  /**
   * Returns the pair orientation suitable for optical duplicates,
   * which always goes by the first then the second end for the strands.
   * This is based on code in MarkDuplicates and ReadEnds.getOrientationByte.
   * Returns one of {@link ReadEnds#RR}, {@link ReadEnds#RF}, {@link ReadEnds#FR}, {@link ReadEnds#FF}
   */
  public byte getOrientationForOpticalDuplicates() {
    final GATKRead read1;
    final GATKRead read2;
    if (first.isFirstOfPair()){
      read1 = first;
      read2 = second;
    } else {
      read1 = second;
      read2 = first;
    }

    final boolean R1R = read1.isReverseStrand();
    final boolean R2R = read2.isReverseStrand();
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
}
