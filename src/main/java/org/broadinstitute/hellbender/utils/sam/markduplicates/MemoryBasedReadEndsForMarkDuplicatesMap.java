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
