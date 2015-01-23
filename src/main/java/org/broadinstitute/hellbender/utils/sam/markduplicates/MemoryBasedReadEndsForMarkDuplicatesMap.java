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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Map from String to ReadEnds object.  Memory-based implementation.  Used for MarkDuplicates.
 *
 * @author alecw@broadinstitute.org
 */
class MemoryBasedReadEndsForMarkDuplicatesMap implements ReadEndsForMarkDuplicatesMap {

    /**
     * Index of this list is sequence index.  Value is map from String {read group id:read name} to ReadEnds.
     * When a ReadEnds is put into this container, it is stored according to the sequenceIndex of the mate
     */
    private List<Map<String, ReadEndsForMarkDuplicates>> mapPerSequence = new ArrayList<Map<String, ReadEndsForMarkDuplicates>>();

    public ReadEndsForMarkDuplicates remove(int mateSequenceIndex, String key) {
        if (mateSequenceIndex >= mapPerSequence.size()) {
            return null;
        }
        return mapPerSequence.get(mateSequenceIndex).remove(key);
    }

    public void put(int mateSequenceIndex, String key, ReadEndsForMarkDuplicates readEnds) {
        while (mateSequenceIndex >= mapPerSequence.size()) {
            mapPerSequence.add(new HashMap<String, ReadEndsForMarkDuplicates>());
        }
        mapPerSequence.get(mateSequenceIndex).put(key, readEnds);
    }

    public int size() {
        int total = 0;
        for (Map<String, ReadEndsForMarkDuplicates> map : mapPerSequence) {
            total += map.size();
        }
        return total;
    }

    /**
     * @return number of elements stored in RAM.  Always <= size()
     */
    public int sizeInRam() {
        return size();
    }
}
