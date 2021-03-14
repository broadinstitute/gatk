package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Map;
import java.util.Set;

final class SplitReadSite {
    private final int position;
    private final Map<String,Integer> sampleCountsMap;

    public SplitReadSite(final int position, final Map<String,Integer> sampleCountsMap) {
        this.position = position;
        this.sampleCountsMap = sampleCountsMap;
    }

    public int getPosition() {
        return position;
    }

    public Map<String,Integer> getSampleCountsMap() {
        return sampleCountsMap;
    }

    public int getSampleCountSum(final Set<String> samples) {
        return sampleCountsMap.entrySet().stream()
                .filter(e -> samples.contains(e.getKey()))
                .mapToInt(e -> e.getValue())
                .sum();
    }

    public double getNormalizedCountSum(final Map<String,Double> sampleCoverageMap) {
        Utils.validateArg(sampleCoverageMap.keySet().containsAll(sampleCountsMap.keySet()), "Coverage missing for one or more samples");
        return sampleCountsMap.entrySet().stream().mapToDouble(e -> e.getValue() / sampleCoverageMap.get(e.getKey())).sum();
    }

    public double getNormalizedCountSum(final Set<String> samples, final Map<String,Double> sampleCoverageMap) {
        Utils.validateArg(sampleCoverageMap.keySet().containsAll(samples), "Coverage missing for one or more samples");
        return sampleCountsMap.entrySet().stream().filter(e -> samples.contains(e.getKey()))
                .mapToDouble(e -> e.getValue() / sampleCoverageMap.get(e.getKey())).sum();
    }

    public int getCountSum() {
        return sampleCountsMap.values().stream().mapToInt(Integer::intValue).sum();
    }
}
