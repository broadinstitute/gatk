package org.broadinstitute.hellbender.utils.read.markduplicates;

/**
 * @author nhomer
 */

/** Stores the minimal information needed for optical duplicate detection. */
public final class PhysicalLocationForMateCigar implements OpticalDuplicateFinder.PhysicalLocation {

    // Information used to detect optical dupes
    short readGroup = -1;
    short tile = -1;
    short x = -1, y = -1;
    short libraryId;

    public PhysicalLocationForMateCigar(final OpticalDuplicateFinder.PhysicalLocation rec) {
        this.setReadGroup(rec.getReadGroup());
        this.setTile(rec.getTile());
        this.setX(rec.getX());
        this.setY(rec.getY());
        this.setLibraryId(rec.getLibraryId());
    }

    @Override
    public short getReadGroup() { return this.readGroup; }

    @Override
    public void setReadGroup(final short rg) { this.readGroup = rg; }

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
    public void setY(final short y) { this.y = y;}

    @Override
    public short getLibraryId() { return this.libraryId; }

    @Override
    public void setLibraryId(final short libraryId) { this.libraryId = libraryId; }

    @Override
    public boolean equals(Object other) {
        if (other instanceof PhysicalLocationForMateCigar) {
            int cmp;
            PhysicalLocationForMateCigar loc = (PhysicalLocationForMateCigar) other;
            cmp = getLibraryId() - loc.getLibraryId();
            if (0 == cmp) cmp = getReadGroup() - loc.getReadGroup();
            if (0 == cmp) cmp = getTile() - loc.getTile();
            if (0 == cmp) cmp = getY() - loc.getY();
            if (0 == cmp) cmp = getX() - loc.getX();
            return 0 == cmp;
        }
        return false;
    }

    @Override
    public int hashCode() {
        int result = getLibraryId();
        result = 31 * result + getReadGroup();
        result = 31 * result + getTile();
        result = 31 * result + getY();
        result = 31 * result + getX();
        return result;
    }
}