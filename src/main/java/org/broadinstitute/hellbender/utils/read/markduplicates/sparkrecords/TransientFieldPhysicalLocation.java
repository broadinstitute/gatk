package org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords;

import picard.sam.util.PhysicalLocation;

/**
 * A common class for holding the fields in PhysicalLocation that we don't want to be serialized by kryo.
 *
 * NOTE: readGroupIndex is not transient as the readgroup is needed in several stages of MarkDuplicatesSpark, but is still
 *       contained in this class to mirror {@link PhysicalLocation}
 */
public abstract class TransientFieldPhysicalLocation extends PairedEnds implements PhysicalLocation {
    // Information used to detect optical dupes
    protected short readGroupIndex = -1;
    protected transient short tile = -1;
    protected transient short x = -1;
    protected transient short y = -1;
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

    // NOTE picard in practice compresses the pixel values to signed shorts for space purposes despite the api using an integer
    @Override
    public void setX(final int x) { this.x = (short)x; }

    @Override
    public int getY() { return this.y; }

    // NOTE picard in practice compresses the pixel values to signed shorts for space purposes despite the api using an integer
    @Override
    public void setY(final int y) { this.y = (short)y; }

    @Override
    public short getLibraryId() { return this.libraryId; }

    @Override
    public void setLibraryId(final short libraryId) { this.libraryId = libraryId; }
}
