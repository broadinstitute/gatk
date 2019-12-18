package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.util.CoordMath;

/**
 * Any class with loci that are potentially on different chromosomes should implement this interface.
 * Consistent with {@link htsjdk.samtools.util.Locatable}, positions should be reported as 1-based and closed at both ends.
 *
 * Note that no particular order is enforced for the two loci.
 */
public interface SVLocatable {

    /**
     * Gets the contig name for first coordinate
     * @return name of the contig
     */
    String getContigA();

    /**
     * @return 1-based first position
     */
    int getPositionA();

    /**
     * Gets the contig name for second coordinate
     * @return name of the contig
     */
    String getContigB();

    /**
     * @return 1-based second position
     */
    int getPositionB();

    /**
     * @return number of bases of reference covered by this interval, or 0 if coordinates on different contigs
     */
    default int getLengthOnReference() {
        if (!isIntrachromosomal()) return 0;
        if (getPositionA() <= getPositionB()) return CoordMath.getLength(getPositionA(), getPositionB());
        return CoordMath.getLength(getPositionB(), getPositionA());
    }

    /**
     * @return true if positions are on the same contig
     */
    default boolean isIntrachromosomal() {
        return getContigA().equals(getContigB());
    }
}
