/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broadinstitute.hellbender.utils.sam.markduplicates;

/**
 * @author nhomer
 */

/** Stores the minimal information needed for optical duplicate detection. */
public class PhysicalLocationForMateCigar implements OpticalDuplicateFinder.PhysicalLocation {

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