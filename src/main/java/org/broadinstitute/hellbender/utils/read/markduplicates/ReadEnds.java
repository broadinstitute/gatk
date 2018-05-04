package org.broadinstitute.hellbender.utils.read.markduplicates;

/** Little struct-like class to hold read pair (and fragment) end data for duplicate marking. */
public abstract class ReadEnds implements picard.sam.util.PhysicalLocation {

    public static final byte F = 0, R = 1, FF = 2, FR = 3, RR = 4, RF = 5;

    public short libraryId;
    public byte orientation;
    public int read1ReferenceIndex = -1;
    public int read1Coordinate = -1;
    public int read2ReferenceIndex = -1;
    public int read2Coordinate = -1;

    // Information used to detect optical dupes
    public short readGroup = -1;
    public short tile = -1;
    public short x = -1, y = -1;

    /** For optical duplicate detection the orientation matters regard to 1st or 2nd end of a mate */
    public byte orientationForOpticalDuplicates = -1;


    public boolean isPaired() { return this.read2ReferenceIndex != -1; }

    @Override
    public short getReadGroup() { return this.readGroup; }

    @Override
    public void setReadGroup(final short readGroup) { this.readGroup = readGroup; }

    @Override
    public short getTile() { return this.tile; }

    @Override
    public void setTile(final short tile) { this.tile = tile; }

    @Override
    public int getX() { return this.x; }

    @Override
    public void setX(final int x) { this.x = (short)x; }

    @Override
    public int getY() { return this.y; }

    @Override
    public void setY(final int y) { this.y = (short)y; }

    @Override
    public short getLibraryId() { return this.libraryId; }

    @Override
    public void setLibraryId(final short libraryId) { this.libraryId = libraryId; }

    /**
     * Returns a single byte that encodes the orientation of the two reads in a pair.
     */
    public static byte getOrientationByte(final boolean read1NegativeStrand, final boolean read2NegativeStrand) {
        if (read1NegativeStrand) {
            if (read2NegativeStrand) return ReadEnds.RR;
            else return ReadEnds.RF;
        } else {
            if (read2NegativeStrand) return ReadEnds.FR;
            else return ReadEnds.FF;
        }
    }
}
