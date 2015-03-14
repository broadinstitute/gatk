package htsjdk.samtools.util;

/**
 * Any class that has a single logical mapping onto the genome should implement Locatable
 * positions should be reported as 1-based and closed at both ends
 *
 */
public interface Locatable {

    /**
     * Gets the contig name for the contig this is mapped to.  May return null if there is no unique mapping.
     * @return name of the contig this is mapped to, potentially null
     */
    String getContig();

    /**
     * @return 1-based start position, undefined if getContig() == null
     */
    int getStart();

    /**
     * @return 1-based closed-ended position, undefined if getContig() == null
     */
    int getEnd();
}
