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
