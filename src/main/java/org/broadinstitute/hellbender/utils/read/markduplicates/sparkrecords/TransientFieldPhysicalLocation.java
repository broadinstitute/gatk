package org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords;

import picard.sam.util.PhysicalLocation;

/**
 * A common interface for holding the fields in PhysicalLocation
 */
public abstract class TransientFieldPhysicalLocation extends PairedEnds implements PhysicalLocation {
    // Information used to detect optical dupes
    protected short readGroupIndex = -1;
    protected transient short tile = -1;
    protected transient int x = -1;
    protected transient int y = -1;
    protected transient short libraryId = -1;

    public TransientFieldPhysicalLocation(int partitionIndex, String name) {
        super(partitionIndex, name);
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
}
