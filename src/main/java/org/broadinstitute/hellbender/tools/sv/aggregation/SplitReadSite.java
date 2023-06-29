package org.broadinstitute.hellbender.tools.sv.aggregation;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Map;
import java.util.Objects;

/**
 * Container class for split read counts for multiple samples at a specific position
 */
public final class SplitReadSite {
    private final String contig;
    private final int position;
    private final boolean strand;
    private final Map<String,Integer> sampleCountsMap;
    private final EvidenceStatUtils.PoissonTestResult result;

    /**
     * @param contig breakpoint contig
     * @param position breakpoint position indicated by the split reads
     * @param strand strandedness of the breakpoint
     * @param sampleCountsMap map with (sample id, split read count > 0) entries
     * @param result statistical test result
     */
    public SplitReadSite(final String contig,
                         final int position,
                         final boolean strand,
                         final Map<String,Integer> sampleCountsMap,
                         final EvidenceStatUtils.PoissonTestResult result) {
        Utils.nonNull(contig);
        Utils.nonNull(sampleCountsMap);
        this.contig = contig;
        this.position = position;
        this.strand = strand;
        this.sampleCountsMap = sampleCountsMap;
        this.result = result;
    }

    public String getContig() {
        return contig;
    }

    public int getPosition() {
        return position;
    }

    public boolean getStrand() {
        return strand;
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
        return position == that.position && strand == that.strand && Objects.equals(sampleCountsMap, that.sampleCountsMap) && Objects.equals(result, that.result);
    }

    @Override
    public int hashCode() {
        return Objects.hash(position, strand, sampleCountsMap, result);
    }
}
