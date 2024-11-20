package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Map;

/**
 * Container class for split read counts for multiple samples at a specific position
 */
final class SplitReadSite {
    private final int position;
    private final Map<String,Integer> sampleCountsMap;

    /**
     * @param position breakpoint position indicated by the split reads
     * @param sampleCountsMap map with (sample id, split read count > 0) entries
     */
    public SplitReadSite(final int position, final Map<String,Integer> sampleCountsMap) {
        Utils.nonNull(sampleCountsMap);
        Utils.validateArg(sampleCountsMap.values().stream().allMatch(c -> c > 0), "Non-positive counts not allowed");
        this.position = position;
        this.sampleCountsMap = sampleCountsMap;
    }

    public int getPosition() {
        return position;
    }

    public boolean hasSample(final String sample) {
        Utils.nonNull(sample);
        return sampleCountsMap.containsKey(sample);
    }

    public int getCount(final String sample) {
        if (hasSample(sample)) {
            return sampleCountsMap.get(sample);
        }
        return 0;
    }
}
