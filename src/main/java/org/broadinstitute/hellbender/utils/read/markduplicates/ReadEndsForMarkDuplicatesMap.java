package org.broadinstitute.hellbender.utils.read.markduplicates;

/**
 * Interface for storing and retrieving ReadEnds objects.  An implementation may be disk-based to
 * reduce memory footprint.
 */
public interface ReadEndsForMarkDuplicatesMap {
    /**
     * Remove element with given key from the map.  Because an implementation may be disk-based,
     * the object returned may not be the same object that was put into the map
     *
     * @param mateSequenceIndex must agree with the value used when the object was put into the map
     * @param key               typically, concatenation of read group ID and read name
     * @return null if the key is not found, otherwise the object removed.
     */
    ReadEndsForMarkDuplicates remove(int mateSequenceIndex, String key);

    /**
     * Store the element in the map with the given key.  It is assumed that the element
     * is not already present in the map.
     *
     * @param mateSequenceIndex use to optimize storage & retrieval.  The same value must be used when trying
     *                          to remove this element.  It is not valid to store the same key with two different mateSequenceIndexes.
     * @param key               typically, concatenation of read group ID and read name
     * @param readEnds          the object to be stored
     */
    void put(int mateSequenceIndex, String key, ReadEndsForMarkDuplicates readEnds);

    /**
     * @return number of elements stored in map
     */
    int size();

    /**
     * @return number of elements stored in RAM.  Always <= size()
     */
    int sizeInRam();
}
