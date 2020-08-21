package org.broadinstitute.hellbender.tools.walkers.sv.GqRecalibrator;

public class MinGq implements Comparable<MinGq> {
    final short minGqHet;
    final short minGqHomVar;

    MinGq(final short minGqHet, final short minGqHomVar) {
        this.minGqHet = minGqHet;
        this.minGqHomVar = minGqHomVar;
    }

    static final MinGq Empty = new MinGq(Short.MIN_VALUE, Short.MIN_VALUE);

    public boolean isEmpty() { return minGqHet == Short.MIN_VALUE && minGqHomVar == Short.MIN_VALUE; }

    @Override
    public String toString() {
        return String.format("(minGqHet:%d, minGqHomVar:%d)", minGqHet, minGqHomVar);
    }

    @Override
    public int compareTo(final MinGq other) {
        // prioritize het, since that's the main intent of the filtering.
        final int hetCompare = Integer.compare(minGqHet, other.minGqHet);
        return hetCompare == 0 ?
            Integer.compare(minGqHomVar, other.minGqHomVar) :
            hetCompare;
    }
}
