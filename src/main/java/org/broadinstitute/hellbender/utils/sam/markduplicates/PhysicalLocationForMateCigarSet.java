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

import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.HashSet;
import java.util.Set;

/**
 * This stores records that are comparable for detecting optical duplicates.
 */
public class PhysicalLocationForMateCigarSet {
    /**
     * We want to return a set of ReadEnds but want to compare based on physical location, hence we store two sets.
     */
    private final Set<ReadEnds> readEnds = new HashSet<ReadEnds>();
    private final Set<PhysicalLocationForMateCigar> physicalLocations = new HashSet<PhysicalLocationForMateCigar>();

    public PhysicalLocationForMateCigarSet() {}

    /** Adds the end to this set, if not already added based on physical location */
    public void add(final ReadEndsForMateCigar end) {
        final PhysicalLocationForMateCigar location = new PhysicalLocationForMateCigar(end);
        if (!physicalLocations.contains(location)) {
            readEnds.add(end);
            physicalLocations.add(new PhysicalLocationForMateCigar(location));
        }
    }

    /** The number of records in this set */
    public int size() { return physicalLocations.size(); }

    /** Removes the end from this set, if present */
    public void remove(final ReadEndsForMateCigar end) {
        final PhysicalLocationForMateCigar location = new PhysicalLocationForMateCigar(end);
        if (physicalLocations.contains(location)) {
            readEnds.remove(end);
            physicalLocations.remove(location);
        }
    }

    /** Gets the set of read ends */
    public Set<ReadEnds> getReadEnds() { return this.readEnds; }

    /** Replaces a given end with the other end.  This ensures that that current is in this set */
    public void replace(final ReadEndsForMateCigar current, final ReadEndsForMateCigar other) {
        final PhysicalLocationForMateCigar location = new PhysicalLocationForMateCigar(current);
        if (!physicalLocations.contains(location)) {
            throw new GATKException("Trying to replace something not in the set");
        }
        this.remove(current);
        this.add(other);
    }
}
