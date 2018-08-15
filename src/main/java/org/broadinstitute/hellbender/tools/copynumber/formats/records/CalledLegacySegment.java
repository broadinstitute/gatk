package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import org.broadinstitute.hellbender.utils.SimpleInterval;

public class CalledLegacySegment extends LegacySegment {

    private CalledCopyRatioSegment.Call call;

    public CalledLegacySegment(final String sampleName, final SimpleInterval interval, final int numProbes, final double segmentMean,
                               final CalledCopyRatioSegment.Call call) {
        super(sampleName, interval, numProbes, segmentMean);
        this.call = call;
    }

    public CalledCopyRatioSegment.Call getCall() {
        return call;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;

        CalledLegacySegment that = (CalledLegacySegment) o;

        return call == that.call;
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + (call != null ? call.hashCode() : 0);
        return result;
    }

    @Override
    public String toString() {
        return "CalledLegacyCopyRatioSegment{" +
                "interval=" + getInterval() +
                ", numPoints=" + getNumProbes() +
                ", meanLog2CopyRatio=" + getSegmentMean() +
                ", call=" + call +
                '}';
    }
}
