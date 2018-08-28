package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

public class CalledLegacySegment extends LegacySegment {

    private final CalledCopyRatioSegment.Call call;

    public CalledLegacySegment(final String sampleName, final SimpleInterval interval, final int numProbes,
                               final double segmentMean,
                               final CalledCopyRatioSegment.Call call) {
        super(sampleName, interval, numProbes, segmentMean);
        Utils.nonNull(call, "Cannot initialize a called legacy segment with a null call.");
        this.call = call;
    }

    public CalledCopyRatioSegment.Call getCall() {
        return call;
    }

    @Override
    public final boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        if (!super.equals(o)) {
            return false;
        }
        final CalledLegacySegment that = (CalledLegacySegment) o;
        return call == that.call;
    }

    @Override
    public final int hashCode() {
        int result = super.hashCode();
        result = 31 * result + call.hashCode();
        return result;
    }

    @Override
    public final String toString() {
        return "CalledLegacySegment{" +
                "interval=" + getInterval() +
                ", numPoints=" + getNumProbes() +
                ", meanLog2CopyRatio=" + getSegmentMean() +
                ", call=" + call +
                '}';
    }
}
