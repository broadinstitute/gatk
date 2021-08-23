package org.broadinstitute.hellbender.tools.sv.aggregation;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Map;
import java.util.Objects;

/**
 * Container class for split read counts for multiple samples at a specific position
 */
public final class SplitReadSite {
    private final int position;
    private final Map<String,Integer> sampleCountsMap;
    private final EvidenceStatUtils.PoissonTestResult result;

    /**
     * @param position breakpoint position indicated by the split reads
     * @param sampleCountsMap map with (sample id, split read count > 0) entries
     */
    public SplitReadSite(final int position, final Map<String,Integer> sampleCountsMap, final EvidenceStatUtils.PoissonTestResult result) {
        Utils.nonNull(sampleCountsMap);
        this.position = position;
        this.sampleCountsMap = sampleCountsMap;
        this.result = result;
    }

    public int getPosition() {
        return position;
    }

    public Double getP() {
        return result == null ? null : result.getP();
    }

    public Double getCarrierSignal() {
        return result == null ? null : result.getCarrierSignal();
    }

    public Double getBackgroundSignal() {
        return result == null ? null : result.getBackgroundSignal();
    }

    public int getCount(final String sample) {
        if (sampleCountsMap.containsKey(sample)) {
            return sampleCountsMap.get(sample);
        }
        return 0;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        SplitReadSite that = (SplitReadSite) o;
        return position == that.position && Objects.equals(sampleCountsMap, that.sampleCountsMap) && Objects.equals(result, that.result);
    }

    @Override
    public int hashCode() {
        return Objects.hash(position, sampleCountsMap, result);
    }
}
