/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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

/** Little struct-like class to hold read pair (and fragment) end data for duplicate marking. */
abstract public class ReadEnds implements OpticalDuplicateFinder.PhysicalLocation {

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
