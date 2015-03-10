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
