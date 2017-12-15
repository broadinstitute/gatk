package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import org.broadinstitute.hellbender.utils.Utils;

public class CalledCopyRatioSegment extends CopyRatioSegment {
    public enum Call {
        AMPLIFICATION("+"),
        DELETION("-"),
        NEUTRAL("0");

        private final String outputString;

        Call(final String outputString) {
            this.outputString = outputString;
        }

        public String getOutputString() {
            return outputString;
        }
    }

    private final Call call;

    public CalledCopyRatioSegment(final CopyRatioSegment segment,
                                  final Call call) {
        super(segment.getInterval(), segment.getNumPoints(), segment.getMeanLog2CopyRatio());
        this.call = Utils.nonNull(call);
    }

    public Call getCall() {
        return call;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        if (!super.equals(o)) {
            return false;
        }
        final CalledCopyRatioSegment that = (CalledCopyRatioSegment) o;
        return call == that.call;
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + call.hashCode();
        return result;
    }

    @Override
    public String toString() {
        return "CalledCopyRatioSegment{" +
                "interval=" + getInterval() +
                ", numPoints=" + getNumPoints() +
                ", meanLog2CopyRatio=" + getMeanLog2CopyRatio() +
                ", call=" + call +
                '}';
    }
}
